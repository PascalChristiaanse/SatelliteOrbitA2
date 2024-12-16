import numpy as np

from Constants import Constants


class PartA:
    def __init__(self):
        self.H_sat = 500  # [km]
        self.e_sat = 0
        self.GPS_transmitter_clock_error = 0.1e-3  # [s]
        self.e_GPS = 0.01
        self.a_GPS = 26560e3  # [m]
        self.R_e = 6473e3  # [m]
        self.mu = 3.986004418e14  # [m^3/s^2]

    def Q1(self):
        """How large are the following effects on the pseudoranges measured by the GPS receiver?
        1) GPS clock offsets (4 points)
        Returns answer in meters
        """
        return Constants.c.value * self.GPS_transmitter_clock_error

    def Q2(self):
        """How large are the following effects on the pseudoranges measured by the GPS receiver?
        2) Light time effect (10 points)
        returns in meters
        """
        a = self.a_GPS
        e = self.e_GPS
        v_gps = np.sqrt(self.mu * (2 / (a * (1 - e)) - (1 / a)))

        # Assume satellites are lined up when signal is sent
        r_sat = self.H_sat + self.R_e
        s_gps_sat = a * (1 - e) - r_sat

        t_sig = s_gps_sat / Constants.c.value
        S_err = t_sig * v_gps
        return S_err

    def Q3(self):
        """How large are the following effects on the pseudoranges measured by the GPS receiver?
        3) Relativistic effect caused by the eccentricity of the GPS orbits (6 points)
        returns in meters
        """
        pass
