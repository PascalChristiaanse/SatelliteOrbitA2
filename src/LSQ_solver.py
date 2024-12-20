import numpy as np
import numba as nb
from typing import Callable

from Constants import Constants




class LSQSolver:
    def __init__(self):
        self.equation = np.array([
            lambda x, u, t: (x[0] - u[0]) / np.sqrt(np.pow(x[0]-u[0], 2) + np.pow(x[1]-u[1], 2) + np.pow(x[2]-u[2], 2)),
            lambda x, u, t: (x[0] - u[0]) / np.sqrt(np.pow(x[0]-u[0], 2) + np.pow(x[1]-u[1], 2) + np.pow(x[2]-u[2], 2)),
            lambda x, u, t: (x[0] - u[0]) / np.sqrt(np.pow(x[0]-u[0], 2) + np.pow(x[1]-u[1], 2) + np.pow(x[2]-u[2], 2)),
            lambda x, u, t: Constants.c.value
        ])

    def construct_jacobian(self, size:int) -> np.ndarray:
        J = np.array([self.equation]*size)
        return J


    def eval_jacobian(self, X: np.ndarray[np.ndarray], U: np.ndarray[np.ndarray]) -> np.ndarray:
        """Evaluate the jacobian matrix for a list of parameters"""
        observ_count = len(X)
        J_evald = np.array([4, observ_count])
        J = self.construct_jacobian(observ_count)
        for i in range(J.shape()[0]):
            for j in range(J.shape()[1]):
                J_evald[i][j] = J[i][j](X[i], U[i], 0)
        return J_evald

    def lsq_iteration():
        pass

def main():
    solver = LSQSolver()

    x_0 = np.array([1,2,3,1,2,3])
    x_1 = np.array([6,5,4,3,2,1])
    X = np.array([x_0, x_1])
    U = np.array([x_1, x_0])
    print(solver.eval_jacobian(X, U))


if __name__ == "__main__":
    main()