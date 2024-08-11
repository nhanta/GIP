# Necessary packages
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from tqdm import tqdm
from utils import (
    normalization, renormalization, 
    rounding, prob, binary_sampler, 
    uniform_sampler, sample_batch_index
)
from sklearn import preprocessing

class GAIN(nn.Module):
    """
    Generator model for GAIN. 
    This neural network takes concatenated data and mask as input and outputs the imputed data.
    """
    def __init__(self, input_dim, h_dim):
        super(GAIN, self).__init__()
        self.fc1 = nn.Linear(input_dim * 2, h_dim)
        self.fc2 = nn.Linear(h_dim, h_dim)
        self.fc3 = nn.Linear(h_dim, input_dim)
    
    def forward(self, x):
        h1 = torch.relu(self.fc1(x))
        h2 = torch.relu(self.fc2(h1))
        return torch.sigmoid(self.fc3(h2))

class Discriminator(nn.Module):
    """
    Discriminator model for GAIN.
    This neural network takes concatenated data and hint as input and outputs the probability of real vs imputed data.
    """
    def __init__(self, input_dim, h_dim):
        super(Discriminator, self).__init__()
        self.fc1 = nn.Linear(input_dim * 2, h_dim)
        self.fc2 = nn.Linear(h_dim, h_dim)
        self.fc3 = nn.Linear(h_dim, input_dim)
    
    def forward(self, x):
        h1 = torch.relu(self.fc1(x))
        h2 = torch.relu(self.fc2(h1))
        return torch.sigmoid(self.fc3(h2))

def round(x):
    """
    Rounding function to convert values to binary (0 or 1).
    """
    return torch.where(x == 0, torch.zeros_like(x), torch.ones_like(x))

def g(X, Y, num_al, N):
    """
    Custom function to compute probability for GAIN based on specific parameters.
    """
    alpha, beta = prob(num_al, N)

    """
    Define a function to compute prob(h_{k+1}|h_1,...,h_k) values for each row pair
    1. Compute intermediate values for the current row of X and all rows of Y
    2. Apply the function to each row of X and Y
    """
    def compute_prob_k(X, y, num_al):
        prob_k = alpha + beta - (alpha + (2 - num_al) * beta) * round(torch.abs(y - X))
        return torch.sum(prob_k, dim=1)
    
    pr = torch.stack([compute_prob_k(X, y, num_al) for y in Y])
    return pr

class Imputation:
    """
    Imputation class to handle missing data imputation using GAIN.

    Args:
        - data_x: Original data with missing values.
        - parameters: Dictionary containing GAIN network parameters:
            - batch_size: Batch size.
            - hint_rate: Hint rate.
            - alpha: Hyperparameter.
            - iterations: Number of iterations.
        - num_al: Number of additional latent variables.
    """
    def __init__(self, data_x, parameters, num_al):
        self.data_x = data_x.T
        self.parameters = parameters
        self.num_al = num_al

    def impute(self, method, num_threads):
        """
        Method to perform data imputation.

        Args:
            - method: String specifying the imputation method ('gip' or 'gain').
            - num_threads: Number of CPU threads to use.
        
        Returns:
            - imputed_data: Data with imputed values.
        """
        torch.set_num_threads(num_threads)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        data_m = 1 - np.isnan(self.data_x)
        batch_size = self.parameters['batch_size']
        hint_rate = self.parameters['hint_rate']
        alpha = self.parameters['alpha']
        iterations = self.parameters['iterations']

        no, dim = self.data_x.shape
        h_dim = dim
        norm_data, norm_parameters = normalization(self.data_x)
        norm_data_x = np.nan_to_num(norm_data, 0)

        # Initialize generator and discriminator models
        generator = GAIN(dim, h_dim).to(device)
        discriminator = Discriminator(dim, h_dim).to(device)

        g_optimizer = optim.Adam(generator.parameters())
        d_optimizer = optim.Adam(discriminator.parameters())

        # Training loop
        for it in tqdm(range(iterations)):
            # Sample batch
            batch_idx = sample_batch_index(no, batch_size)
            X_mb = norm_data_x[batch_idx, :]
            M_mb = data_m[batch_idx, :]
            # Sample random vectors  
            Z_mb = uniform_sampler(0, 0.01, batch_size, dim)
            Z_mb = preprocessing.normalize(Z_mb)

            # Sample hint vectors
            H_mb_temp = binary_sampler(hint_rate, batch_size, dim)
            H_mb = M_mb * H_mb_temp

            # Combine random vectors with observed vectors
            X_mb = M_mb * X_mb + (1 - M_mb) * Z_mb
            X_mb, M_mb, H_mb = map(lambda x: torch.tensor(x, dtype=torch.float32).to(device), [X_mb, M_mb, H_mb])

            if method == 'gip':
                G_sample = generator(torch.cat([X_mb, M_mb], dim=1))
            elif method == 'gain':
                G_sample = generator(torch.cat([X_mb, M_mb], dim=1))
            else:
                raise ValueError("Method is wrong, please choose a correct method")

            Hat_X = X_mb * M_mb + G_sample * (1 - M_mb)
            D_prob = discriminator(torch.cat([Hat_X, H_mb], dim=1))

            D_loss = -torch.mean(M_mb * torch.log(D_prob + 1e-8) + (1 - M_mb) * torch.log(1. - D_prob + 1e-8))
            G_loss = -torch.mean((1 - M_mb) * torch.log(D_prob + 1e-8)) + alpha * torch.mean((M_mb * X_mb - M_mb * G_sample) ** 2) / torch.mean(M_mb)

            d_optimizer.zero_grad()
            D_loss.backward()
            d_optimizer.step()

            g_optimizer.zero_grad()
            G_loss.backward()
            g_optimizer.step()

        # Impute missing data
        Z_mb = uniform_sampler(0, 0.01, no, dim)
        Z_mb = preprocessing.normalize(Z_mb)
        X_mb = norm_data_x
        X_mb = data_m * X_mb + (1 - data_m) * Z_mb
        X_mb, M_mb = map(lambda x: torch.tensor(x, dtype=torch.float32).to(device), [X_mb, data_m])

        with torch.no_grad():
            if method == 'gip':
                imputed_data = generator(torch.cat([X_mb, M_mb], dim=1)).cpu().numpy()
            elif method == 'gain':
                imputed_data = generator(torch.cat([X_mb, M_mb], dim=1)).cpu().numpy()

        imputed_data = data_m * norm_data_x + (1 - data_m) * imputed_data
        imputed_data = renormalization(imputed_data, norm_parameters)
        imputed_data = rounding(imputed_data.T, self.data_x.T)

        return imputed_data
