# Necessary packages
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from tqdm import tqdm
from sklearn import preprocessing

from utils import rounding, normalization, renormalization
from utils import prob, binary_sampler, uniform_sampler, sample_batch_index

# Define the activation function G
def G(X, Y, num_alleles, num_haplotypes):

    alpha, beta = prob(num_alleles, num_haplotypes)

    # Helper function to round tensor values to 0 or 1
    def round(x):
        return torch.where(x == 0, torch.zeros_like(x), torch.ones_like(x))

    # Compute prob(h_{k+1,j}|h_1,...,h_k) for each row pair.
    def compute_prob_k(X, y):
        # prob_k = alpha + beta - alpha * round(torch.abs(y - X))
        gamma = 10 ** 4
        #prob_k = alpha * torch.exp(-gamma * (y - X) ** 2) + beta
        prob_k =  alpha / (1 + gamma * (y - X) ** 2) + beta 
        return (1/num_haplotypes) * torch.sum(prob_k, dim=0) 

    pr = torch.stack([compute_prob_k(X, y) for y in Y], dim=0)
    return pr

# Discriminator Network
class Discriminator(nn.Module):
    def __init__(self, num_features):
        super(Discriminator, self).__init__()
        self.fc1 = nn.Linear(2 * num_features, num_features)
        self.fc2 = nn.Linear(num_features, num_features)
        self.fc3 = nn.Linear(num_features, num_features)

    def forward(self, x, h):
        inputs = torch.cat([x, h], dim=1)
        h1 = torch.relu(self.fc1(inputs))
        h2 = torch.relu(self.fc2(h1))
        logits = self.fc3(h2)
        return torch.sigmoid(logits)

# GIP Generator Network
class GIP_Generator(nn.Module):
    def __init__(self, num_features, num_alleles):
        super(GIP_Generator, self).__init__()
        self.fc1 = nn.Linear(2 * num_features, num_features)
        self.fc2 = nn.Linear(num_features, num_features)
        self.fc3 = nn.Linear(num_features, num_features)
        self.num_alleles = num_alleles

    def forward(self, x, m):
        inputs = torch.cat([x, m], dim=1)
        z1 = self.fc1(inputs)
        h1 = torch.relu(z1)
        z2 = self.fc2(h1)
        #h2 = torch.relu(z2)
        h2 = nn.functional.gumbel_softmax(z2, tau = 1.0, hard = True)
        z3 = self.fc3(h2)
        h3 = G(h2, z3, self.num_alleles, x.shape[0])
        #h2 = G(h1, z1, self.num_alleles, x.shape[0])
        #logits = self.fc3(h2)
        return h3
        #torch.sigmoid(logits)

# GAIN Generator Network
class GAIN_Generator(nn.Module):
    def __init__(self, num_features):
        super(GAIN_Generator, self).__init__()
        self.fc1 = nn.Linear(2 * num_features, num_features)
        self.fc2 = nn.Linear(num_features, num_features)
        self.fc3 = nn.Linear(num_features, num_features)

    def forward(self, x, m):
        inputs = torch.cat([x, m], dim=1)
        h1 = torch.relu(self.fc1(inputs))
        h2 = torch.relu(self.fc2(h1))
        logits = self.fc3(h2)
        return torch.sigmoid(logits)

# GIP Class - Main Class for Missing Data Imputation
class GIP:
    '''Impute missing values in data_x
    
    Args:
      - data_x: original data with missing values
      - parameters: Imputation parameters, including:
        - batch_size: Batch size for mini-batch training
        - hint_rate: Rate for adding hints to discriminator
        - alpha: Hyperparameter to balance the generator loss
        - iterations: Number of training iterations
      - num_alleles: Number of alleles 
    
    Returns:
      - imputed_data: Data matrix with missing values imputed
    '''
    def __init__(self, data_x, gip_parameters, system_parameters):
        self.data_x = data_x.T
        self.gip_parameters = gip_parameters
        self.system_parameters = system_parameters  
    
    def impute(self, method):
        ncpu = self.system_parameters['num_cpus']
        ngpu = self.system_parameters['num_gpus']
        device = torch.device('cuda' if ngpu > 0 and torch.cuda.is_available() else 'cpu')
        # Set the number of threads for parallelism in PyTorch
        if ncpu: torch.set_num_threads(ncpu)
        # Check for available GPUs and set device
        print(f"Using device: {device}")

        # Mask matrix indicating observed (1) and missing (0) values
        data_m = 1 - np.isnan(self.data_x)

        # Extract system parameters
        batch_size = self.gip_parameters['batch_size']
        hint_rate = self.gip_parameters['hint_rate']
        alpha = self.gip_parameters['alpha']
        iterations = self.gip_parameters['iterations']
        num_alleles = self.gip_parameters['num_alleles']
        # Get the number of samples (N) and features (D)
        num_haplotypes, num_features = self.data_x.shape

        # Normalize the data
        parameters = {'min_val': 0, 'max_val': num_alleles - 1}
        norm_data, norm_parameters = normalization(self.data_x, parameters)
        norm_data_x = np.nan_to_num(norm_data, 0)

        # Convert data and mask to PyTorch tensors and move to device
        norm_data_x = torch.tensor(norm_data_x, dtype=torch.float32).to(device)
        data_m = torch.tensor(data_m, dtype=torch.float32).to(device)

        # Initialize Discriminator and appropriate Generator model
        discriminator = Discriminator(num_features).to(device)
        if method == "gip":
            generator = GIP_Generator(num_features, num_alleles).to(device)
        elif method == "gain":
            generator = GAIN_Generator(num_features).to(device)
        else:
            raise ValueError("Invalid method. Choose 'gip_edit' or 'gain'.")

        # If using multiple GPUs, wrap the models in DataParallel
        if ngpu > 1:
            discriminator = nn.DataParallel(discriminator, device_ids=range(ngpu))
            generator = nn.DataParallel(generator, device_ids=range(ngpu))

        # Optimizers for both models
        d_optimizer = optim.Adam(discriminator.parameters())
        g_optimizer = optim.Adam(generator.parameters())

        # Loss function - Binary Cross-Entropy (BCE) for classification
        #bce_loss = nn.BCELoss(reduction='mean')

        # Training loop
        for it in tqdm(range(iterations)):    
            # Sample a batch of data
            batch_idx = sample_batch_index(num_haplotypes, batch_size)
            X_mb = norm_data_x[batch_idx, :]  # Batch of normalized data
            M_mb = data_m[batch_idx, :]       # Corresponding mask

            # Sample random noise for missing values
            Z_mb = uniform_sampler(0, 0.01, batch_size, num_features)
            Z_mb = torch.tensor(preprocessing.normalize(Z_mb), dtype=torch.float32).to(device)

            # Sample hint vectors (to help discriminator)
            H_mb_temp = binary_sampler(hint_rate, batch_size, num_features)
            H_mb = M_mb * torch.tensor(H_mb_temp, dtype=torch.float32).to(device)
            
            # Combine random vectors with observed data
            X_mb = M_mb * X_mb + (1 - M_mb) * Z_mb 

            # Discriminator optimization step
            G_sample = generator(X_mb, M_mb)  # Generate data using the generator
            Hat_X = X_mb * M_mb + G_sample * (1 - M_mb)  # Combine generated and observed data
            D_prob = discriminator(Hat_X, H_mb)  # Pass through discriminator

            # Compute discriminator loss
            D_loss = -torch.mean(M_mb * torch.log(D_prob + 1e-8) + (1 - M_mb) * torch.log(1. - D_prob + 1e-8))
            d_optimizer.zero_grad()
            D_loss.backward()  # Backpropagation
            d_optimizer.step()

            # Generator optimization step
            G_sample = generator(X_mb, M_mb)  # Generate data
            Hat_X = X_mb * M_mb + G_sample * (1 - M_mb)  # Combine with observed data
            D_prob = discriminator(Hat_X, H_mb)  # Discriminator output

            # Generator loss
            G_loss_temp = -torch.mean((1 - M_mb) * torch.log(D_prob + 1e-8))
            MSE_loss = torch.mean((M_mb * X_mb - M_mb * G_sample) ** 2) / torch.mean(M_mb)
            G_loss = G_loss_temp + alpha * MSE_loss  # Combine losses

            g_optimizer.zero_grad()
            G_loss.backward()  # Backpropagation
            g_optimizer.step()

        # Final imputation using trained generator
        Z_mb = uniform_sampler(0, 0.01, num_haplotypes, num_features)
        X_mb = norm_data_x
        Z_mb = torch.tensor(preprocessing.normalize(Z_mb), dtype=torch.float32).to(device)
        X_mb = data_m * X_mb + (1 - data_m) * Z_mb 

        with torch.no_grad():
            G_sample = generator(X_mb, data_m)  # Generate final sample

        imputed_data = data_m * norm_data_x + (1 - data_m) * G_sample  # Imputed data
        imputed_data = imputed_data.cpu().numpy()  # Move to CPU for conversion to numpy

        # Denormalize and round imputed data
        imputed_data = renormalization(imputed_data, norm_parameters)  
        imputed_data = rounding(imputed_data.T, self.data_x.T)  # Apply rounding to match original format
        
        return imputed_data
    
    