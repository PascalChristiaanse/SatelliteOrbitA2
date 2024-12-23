import numpy as np
from numpy.linalg import inv, norm
from numpy import cos, sin

from Constants import Constants
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))
from data.read_txt_python import Data

from OrbitPlotter import OrbitPlotter


class LSQSolver:
    def __init__(self, equations, derivatives):
        self.equations = equations
        self.derivatives = derivatives

    @staticmethod
    def get_covariance_matrix(std, corr, n):
        matrix = corr * np.ones((n, n))

        # Set the diagonal to 1 (perfect correlation with itself)
        np.fill_diagonal(matrix, 1)
        matrix = matrix * np.pow(std, 2)
        return matrix

    def construct_information(self, size: int):
        A = np.array([self.equations] * size)
        return A

    def construct_jacobian(self, size: int) -> np.ndarray:
        J = np.array([self.derivatives] * size)
        return J

    def eval_information(
        self, x_bar: np.ndarray, U: np.ndarray[np.ndarray]
    ) -> np.ndarray[np.ndarray]:
        """Evaluate the jacobian matrix for a list of parameters

        Args:
            X (np.ndarray): estimate of parameters
            U (np.ndarray[np.ndarray]): equation parameters

        Returns:
            np.ndarray: Jacobian evaluated at X for U
        """

        A = self.construct_information(U.shape[0])
        A_evald = np.zeros(A.shape)
        for i in range(A.shape[0]):
            A_evald[i] = A[i](x_bar, U[i], 0)
        return A_evald

    def eval_jacobian(
        self, x_bar: np.ndarray, U: np.ndarray[np.ndarray]
    ) -> np.ndarray[np.ndarray]:
        """Evaluate the jacobian matrix for a list of parameters

        Args:
            X (np.ndarray): estimate of parameters
            U (np.ndarray[np.ndarray]): equation parameters

        Returns:
            np.ndarray: Jacobian evaluated at X for U
        """

        J = self.construct_jacobian(U.shape[0])
        J_evald = np.zeros(J.shape)
        for i in range(J.shape[0]):
            for j in range(J.shape[1]):
                J_evald[i][j] = J[i][j](x_bar, U[i], 0)
        return J_evald

    def lsq_iteration(
        self,
        x_bar: np.ndarray,
        y_bar: np.ndarray,
        U: np.ndarray[np.ndarray],
        P_yy: np.ndarray[np.ndarray],
    ):
        """Performs a single non-linear LSQ estimate

        Args:
            x_bar (np.ndarray): initial guess of parameters required to match observations
            y_bar (np.ndarray): observations
            U (np.ndarray[np.ndarray]): Parameter matrix for jacobian
            covariance (np.ndarray[np.ndarray]): covariance matrix of observations

        Returns:
            _type_: updated estimate
        """
        J = self.eval_jacobian(x_bar=x_bar, U=U)  # A matrix at x_bar
        f_x = self.eval_information(x_bar=x_bar, U=U)
        # print(J)
        # raise KeyError
        estimate = x_bar + inv(J.T @ inv(P_yy) @ J) @ J.T @ inv(P_yy) @ (y_bar - f_x)
        return estimate

    def solve(
        self,
        guess: np.ndarray,
        observation: np.ndarray,
        parameters: np.ndarray[np.ndarray],
        std: float,
        cofac: float,
    ):
        estimate = guess
        P_yy = self.get_covariance_matrix(std, cofac, parameters.shape[0])
        for i in range(1000):
            new_estimate = self.lsq_iteration(
                x_bar=estimate, y_bar=observation, U=parameters, P_yy=P_yy
            )
            if np.linalg.norm(estimate - new_estimate) < 1e-6:
                # print("converged at:", new_estimate, "after", i, "Iterations")
                solution = new_estimate
                return solution
            else:
                estimate = new_estimate
        print("DIDNT converge")


def main():
    data = Data()

    # Get reference orbit
    ref_sol = []
    for i in range(len(data.rx)):
        ref_sol.append(
            [
                data.rx[i] + 0.00707 * data.vx[i],
                data.ry[i] + 0.00707 * data.vy[i],
                data.rz[i] + 0.00707 * data.vz[i],
                data.vx[i],
                data.vy[i],
                data.vz[i],
            ]
        )
    ref_sol = np.array(ref_sol)

    # Calculate LSQ solution
    # Part D

    # fmt: off
    rho = lambda x, u, t: norm(x[0:3]-u[0:3]) + (Constants.c.value/1000)*(x[3]-u[3])  # noqa: E731
    derivatives = np.array(
        [
            lambda x, u, t: (x[0] - u[0])/ norm(x[0:3]-u[0:3]),
            lambda x, u, t: (x[1] - u[1])/ norm(x[0:3]-u[0:3]),
            lambda x, u, t: (x[2] - u[2])/ norm(x[0:3]-u[0:3]),
            lambda x, u, t: Constants.c.value/1000,
        ]
    )

    # fmt: on
    solver = LSQSolver(equations=rho, derivatives=derivatives)

    rx0 = data.rx_gps[0] - data.clk_gps[0] * data.vx_gps[0]
    ry0 = data.ry_gps[0] - data.clk_gps[0] * data.vy_gps[0]
    rz0 = data.rz_gps[0] - data.clk_gps[0] * data.vz_gps[0]
    clk0 = data.clk_gps[0]

    # ca_range0 = data.CA_range[0]
    U0 = np.column_stack((rx0, ry0, rz0, clk0))
    U0 = U0[~np.all(U0 == 0, axis=1)]
    # y_0 = ca_range0[ca_range0 != 0]
    x_0 = np.array([6374, 0, 0, 0.001])
    derivatives[0](x_0, U0[0], 0)

    # Part E
    lsq_solution = []
    for i in range(200):
        rx = data.rx_gps[i] - data.clk_gps[i] * data.vx_gps[i]
        ry = data.ry_gps[i] - data.clk_gps[i] * data.vy_gps[i]
        rz = data.rz_gps[i] - data.clk_gps[i] * data.vz_gps[i]
        clk = data.clk_gps[i]
        ca_range = data.CA_range[i]
        U = np.column_stack((rx, ry, rz, clk))
        U = U[~np.all(U == 0, axis=1)]
        y_i = ca_range[ca_range != 0]

        x_guess = np.array([6374, 0, 0, 0.001])

        solution = solver.solve(x_guess, y_i, U, 3 / 1000, 0.2)

        lsq_solution.append([*solution, data.t[i] + 0.00707])
    lsq_solution = np.array(lsq_solution)

    # Plotter stuff
    plotter = OrbitPlotter()
    plotter.add_orbit(ref_sol, name="reference", color="blue")
    plotter.add_orbit(lsq_solution[:, :-1], name="LSQ Solution", color="green")
    plotter.plot()

    from ResidualPlotter import ResidualPlotter

    residual = ResidualPlotter()
    residual.add_line_plot(
        data.t - 959300000,
        1000 * np.linalg.norm(ref_sol[:, :-3] - lsq_solution[:, :-2], axis=1),
        name="Residual lsq vs reference",
        color="red",
    )
    # residual.plot()

    print(["{0:0.2f}".format(i * 1000) for i in lsq_solution[0]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_solution[1]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_solution[2]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_solution[3]])

    # Part F
    R_z = lambda g: np.array(  # noqa: E731
        [
            [cos(g), sin(g), 0],
            [-sin(g), cos(g), 0],
            [0, 0, 1],
        ]
    )
    r_t = lambda r, v, tau: R_z(Constants.omg_e.value * tau) @ (r - tau * v)  # noqa: E731
    relv = lambda r, v: (-2 / (Constants.c.value / 1000)) * np.dot(r, v)  # noqa: E731
    rho_lt = (  # noqa: E731
        lambda x, u, t: norm(x[0:3] - r_t(u[0:3], u[3:6], u[7]))
        + (Constants.c.value / 1000) * (x[3] - u[6])
        - relv(u[0:3], u[3:6])
    )

    # fmt: off
    jaco_lt = np.array(
        [
            lambda x, u, t: (rho_lt(x+np.array([1e-6,0,0,0]),u,t)-rho_lt(x-np.array([1e-6,0,0,0]),u,t))/2e-6,
            lambda x, u, t: (rho_lt(x+np.array([0,1e-6,0,0]),u,t)-rho_lt(x-np.array([0,1e-6,0,0]),u,t))/2e-6,
            lambda x, u, t: (rho_lt(x+np.array([0,0,1e-6,0]),u,t)-rho_lt(x-np.array([0,0,1e-6,0]),u,t))/2e-6,
            lambda x, u, t: Constants.c.value/1000,
        ]
    )

    # fmt: on

    solver_lt_rel = LSQSolver(equations=rho_lt, derivatives=jaco_lt)
    lsq_lt_sol = []
    for i in range(200):
        rx = data.rx_gps[i] - data.clk_gps[i] * data.vx_gps[i]
        ry = data.ry_gps[i] - data.clk_gps[i] * data.vy_gps[i]
        rz = data.rz_gps[i] - data.clk_gps[i] * data.vz_gps[i]
        vx = data.vx_gps[i]
        vy = data.vy_gps[i]
        vz = data.vz_gps[i]
        clk = data.clk_gps[i]
        ca_range = data.CA_range[i]
        y_i = ca_range[ca_range != 0]

        U = np.column_stack(
            (rx, ry, rz, vx, vy, vz, clk, ca_range / (Constants.c.value / 1000))
        )
        U = U[~np.all(U == 0, axis=1)]
        x_guess = np.array([6374, 0, 0, 0.001])
        solution = solver_lt_rel.solve(x_guess, y_i, U, 3 / 1000, 0.2)
        lsq_lt_sol.append([*solution, data.t[i] + 0.00707])

    lsq_lt_sol = np.array(lsq_lt_sol)
    residual.add_line_plot(
        data.t - 959300000,
        1000 * np.linalg.norm(ref_sol[:, :-3] - lsq_lt_sol[:, :-2], axis=1),
        name="Residual LightT+Relativity vs Ref",
        color="green",
    )
    residual.plot()
    from TimeOffsetPlotter import TimePlotter

    timeplotter = TimePlotter()
    timeplotter.add_line_plot(
        data.t - 959300000,
        lsq_lt_sol[:, 3]*1000,
        name="Receiver Clockoffset lt+rel",
        color="black",
    )
    timeplotter.plot()
    print("Relativistic and lighttime effects")
    print(["{0:0.2f}".format(i * 1000) for i in lsq_lt_sol[0]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_lt_sol[1]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_lt_sol[2]])
    print(["{0:0.2f}".format(i * 1000) for i in lsq_lt_sol[3]])

    std_sol = np.sqrt(np.sum(norm(lsq_lt_sol[:, :-2]-ref_sol[:, :-3], axis=1))/(len(lsq_lt_sol)-1))
    print("standard deviation:", std_sol*1000)
    print("PDOP:", np.sqrt(3*pow(std_sol*1000,2)))
if __name__ == "__main__":
    main()
