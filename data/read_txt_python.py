import numpy as np
import pandas as pd
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))



class Data:
    def __init__(self):
        self.t = np.loadtxt('data/t.txt') # epochs - s

        self.CA_range = np.loadtxt('data/CA_range.txt') # pseudorange observations from CA code - km

        self.PRN_ID = np.loadtxt('data/PRN_ID.txt') # PRN ID of tracked GPS satellites

        self.clk_gps = np.loadtxt('data/clk_gps.txt') # clock correction for GPS sats (transmitters) - s

        self.rx_gps = np.loadtxt('data/rx_gps.txt') # GPS satellite positions (transmitters) - km
        self.ry_gps = np.loadtxt('data/ry_gps.txt')
        self.rz_gps = np.loadtxt('data/rz_gps.txt')

        self.vx_gps = np.loadtxt('data/vx_gps.txt') # GPS satellite velocities (transmitters) - km/s
        self.vy_gps = np.loadtxt('data/vy_gps.txt')
        self.vz_gps = np.loadtxt('data/vz_gps.txt')

        self.rx = np.loadtxt('data/rx.txt') # precise positions (receivers) - km
        self.ry = np.loadtxt('data/ry.txt')
        self.rz = np.loadtxt('data/rz.txt')

        self.vx = np.loadtxt('data/vx.txt') # precise velocities (receivers) - km/s
        self.vy = np.loadtxt('data/vy.txt')
        self.vz = np.loadtxt('data/vz.txt')
   
    
    
def main():
    data = Data()  
        
if __name__ == "__main__":
    main()