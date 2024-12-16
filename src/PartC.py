import numpy as np

from Constants import Constants


class PartC:
    def __init__(self):
        self.std = 3  # [m]
        self.correlation = 0.2

    def get_covariance_matrix(self, std, corr, n):
        matrix = corr * np.ones((n, n))

        # Set the diagonal to 1 (perfect correlation with itself)
        np.fill_diagonal(matrix, 1)
        matrix = matrix * np.pow(std, 2)
        return matrix

    def Q1(self):
        return """The covariance matrix is an nxn matrix with every cell equal to the correlation factor,
        except for the main diagonal on which all values are equal to 1. 
        This base matrix is multiplied with the standard deviation squared (also called the variance)
        Please note that this procedure is exclusively valid if all correlations are the same, and all standard deviations are the same.
        We find the matrix for std=3, correlation factor=0.2, n=4 (4 observations) \n""" + str(self.get_covariance_matrix(self.std, self.correlation, 4))
