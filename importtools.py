# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 14:13:47 2025

@author: aengstrom
"""
import netCDF4 as ncdf
import numpy as np

class Chromatogram:
    def __init__(self, filename, dataformat):
        self.format = dataformat
        self.filename = filename
        self.chromatogram = self.generatechrom()
        self.peaks = self.generateareas()
        


    def readcdf(self):
        # Open the CDF file
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            print("Global attributes:", rootgrp.ncattrs())
            
            # Extract all variables
            print("Variables in file:", rootgrp.variables.keys())
            # Get signal (ordinate values)
            print(rootgrp.variables["actual_sampling_interval"])

    def generatechrom(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            if "ordinate_values" in rootgrp.variables:
                signal = np.array(rootgrp.variables["ordinate_values"][:])
            else:
                raise KeyError("ordinate_values not found in the file!")
            if "actual_run_time_length" in rootgrp.variables and "actual_delay_time" in rootgrp.variables:
                runtime, starttime = rootgrp.variables["actual_run_time_length"][:], rootgrp.variables["actual_delay_time"]
                rt = np.linspace(start=starttime, stop=runtime, num=len(signal))
            else:
               raise KeyError("runtime values not found in the file!")
            return np.vstack((rt, signal))
                
    def generateareas(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            peak_names = [
                b"".join(np.ma.filled(row, b"")).decode("utf-8").strip()
                for row in rootgrp.variables["peak_name"][:]
            ]
            peak_amount = rootgrp.variables["peak_amount"][:]
            return np.array((peak_names,peak_amount))
            
        

if __name__ == "__main__":
    c1 = Chromatogram("C:/Users/aengstrom/Downloads/rbch06cdat-Front Signal.cdf", "m")


