# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 14:13:47 2025

@author: aengstrom
"""
#%%
import logging
import netCDF4 as ncdf
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
import re
from collections import defaultdict
from dateutil import parser
import os
from chemformula import ChemFormula
import pubchempy as pcp

class Chromatogram:
    """Handles chromatogram data from .cdf files."""
    # bp_voc_list = [ChemFormula('CH3CH2CH2CH2CH2CH3', name = 'n-hexane', cas = 110_54_3),
    #                ChemFormula('C5H9(CH3)', name = 'methylcyclopentane', cas = 96_37_7),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),
    #                ChemFormula(),]
    def __init__(self, filename, dataformat="cdf"):
        self.format = dataformat
        self.filename = Path(filename)
        self._datetime = None
        self._chromatogram = None
        self._peakamounts = None
        self._peakwindows = None
        self._peaklocations = None
        self._loaded = False
        self._error = None

    def load(self):
        """Explicitly load all data"""
        try:
            self._datetime = self._get_datetime()
            self._chromatogram = self._generate_chrom()
            self._peakamounts = self._generate_class_attributes('peakamounts')
            self._peakwindows = self._generate_class_attributes('peakwindows')
            self._peaklocations = self._generate_class_attributes('peaklocations')
            self._loaded = True
            self._error = None
        except Exception as e:
            self._error = str(e)
            raise

    @property
    def datetime(self):
        if not self._loaded and self._datetime is None:
            self._datetime = self._get_datetime()
        return self._datetime

    @property
    def chromatogram(self):
        if not self._loaded and self._chromatogram is None:
            self._chromatogram = self._generate_chrom()
        return self._chromatogram

    @property
    def peakamounts(self):
        if not self._loaded and self._peakamounts is None:
            self._peakamounts = self._generate_class_attributes('peakamounts')
        return self._peakamounts
    
    @property
    def peakwindows(self):
        if not self._loaded and self._peakwindows is None:
            self._peakwindows = self._generate_class_attributes('peakwindows')
        return self._peakwindows

    @property
    def peaklocations(self):
        if not self._loaded and self._peaklocations is None:
            self._peaklocations = self._generate_class_attributes('peaklocations')
        return self._peaklocations

        
    def _get_datetime(self):
        """Private method to get datetime"""
        try:
            with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
                # Example: get datetime from file metadata
                if 'dataset_date_time_stamp' in rootgrp.ncattrs():
                    date_string = rootgrp.dataset_date_time_stamp
                    return parser.parse(date_string)
                else:
                    # Fallback: use file modification time
                    return datetime.fromtimestamp(self.filename.stat().st_mtime)
        except Exception as e:
            print(f"Error getting datetime: {e}")
            return None
        
    def _generate_chrom(self):
        """Private method to generate chromatogram"""
        try:
            with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
                required_vars = ["ordinate_values", "actual_run_time_length", "actual_delay_time", 'actual_sampling_interval']
                for var in required_vars:
                    if var not in rootgrp.variables:
                        raise KeyError(f"Missing '{var}' in {self.filename}")
                
                signal = np.array(rootgrp.variables["ordinate_values"][:])
                runtime = float(rootgrp.variables["actual_run_time_length"][0])
                starttime = float(rootgrp.variables["actual_delay_time"][0])
                interval = float(rootgrp.variables['actual_sampling_interval'][0])
                rt = np.arange(starttime, runtime, interval)
                return np.vstack((rt, signal))
        except Exception as e:
            print(f"Error accessing file {self.filename}: {e}")
            return None
            
    def _generate_class_attributes(self, attribute: str):
        """Private method to generate attributes"""
        attribute_guide = {'peaklocations': ["peak_name", 'baseline_start_time', 'baseline_stop_time', 'baseline_start_value', 'baseline_stop_value', "peak_retention_time"],
                           'peakamounts': ["peak_name", "peak_amount"],
                           'peakwindows': ["peak_name", "peak_start_time", "peak_end_time", "peak_retention_time"]}
        try:
            with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
                required_vars = attribute_guide[attribute]
                for var in required_vars:
                    if var not in rootgrp.variables:
                        raise KeyError(f"Missing '{var}' in {self.filename}")
                data = {var: ncdf.chartostring(rootgrp.variables[var][:]) if var == "peak_name" 
                        else rootgrp.variables[var][:]
                        for var in required_vars}

                return pd.DataFrame(data)
        except Exception as e:
            print(f"Error reading peak start or end times: {e}")
    
    def list_netcdf_variables(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return rootgrp.variables.keys()
    def list_netcdf_attributes(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return rootgrp.__dict__
    def examine_netcdf_variable(self, variable):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return rootgrp.variables[variable][:], rootgrp.variables[variable].__dict__
    def examine_netcdf_attribute(self, attribute):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return getattr(rootgrp, attribute)
    

    def __repr__(self):
        return f"Chromatogram({self.filename.name})"
 
class DataAnalysis:
    """Represents a single year, month, or day of chromatogram data."""

    def __init__(self, site_path, year, month = None, day = None):
        """
        Creates DataAnalysis object that contains chromatogram objects collected over a year, month, or day.
        Parameters:
            site_path: str | Path -> path to site root
            year: str or int-> year of interest. If no other date inputs are entered, data from the whole year will be extracted
            month: str or int -> month of interest.
            day: str or int -> day of interest
        """
        self.date = [year, month, day]
        self.folder = self._find_root_dir(site_path = site_path, date = self.date)
        self.filename = self.folder.stem
        self.samples = self._get_chromatograms(self.folder, 's')
        self.cvs = self._get_chromatograms(self.folder, 'c')
        self.lcs = self._get_chromatograms(self.folder, 'e')
        self.blank = self._get_chromatograms(self.folder, 'b')
        self.rts = self._get_chromatograms(self.folder, 'q')
        
    def _find_root_dir(self, site_path, date: list):
        root_dir = Path(site_path) / str(date[0])
        for j in date[1:]:
            if j:
                root_dir = root_dir / f"{int(j):02d}"
        return root_dir
                
            
    def _chrom_type(self, filename):
        pattern = r"(?P<site>[A-Z]{2})(?P<sample_type>[A-Z])(?P<month>[A-Z])(?P<day>\d{2})(?P<hour>[A-Z]).*-(?P<column>Front Signal|Back Signal)"
        match = re.match(pattern, filename, re.IGNORECASE)
        if match:
            return match.groupdict()
        else:
            print(f"Could not identify match for {filename}")
            return {
                "site": None,
                "sample_type": None,
                "month": None,
                "day": None,
                "hour": None,
                "column": None
            }
     
    def _get_chromatograms(self, path, sample_type_letter):
        """
        Returns a list of Sample/Blank/RTS/CVS/LCS objects
        by pairing Front and Back Signal files.
        """

        if not path.exists():
            print(f"Warning: Folder does not exist → {path}")
            return []

        # Step 1: scan all files and parse metadata
        files_by_run = defaultdict(dict)
        for file in path.rglob("**/*.cdf"):
            info = self._chrom_type(file.stem)
            if not info["sample_type"]:
                continue  # skip unrecognized files

            # unique key for a run: site+sample_type+month+day+hour
            run_key = (info["site"], info["sample_type"], info["month"], info["day"], info["hour"])
            
            if info["column"].lower() == "front signal":
                files_by_run[run_key]["front"] = file
            elif info["column"].lower() == "back signal":
                files_by_run[run_key]["back"] = file

        # Step 2: instantiate objects based on sample_type
        result_objects = []
        type_map = {
            "s": Sample,
            "b": Blank,
            "c": CVS,
            "q": RTS,
            "e": LCS,
            "d": DetectionLimit,
            "m": Calibration
        }

        for run_key, paths in files_by_run.items():
            if run_key[1] == sample_type_letter:
                front = paths.get("front")
                back = paths.get("back")
                if front and back:
                    sample_type = run_key[1].lower()
                    cls = type_map.get(sample_type)
                    if cls:
                        result_objects.append(cls(front, back))
                    else:
                        print(f"⚠️ Unknown sample type '{sample_type}' for run {run_key}")
                else:
                    print(f"⚠️ Missing front or back file for run {run_key}")

        return result_objects
        
        
    def __str__(self):
        return f"{self.folder}"
        

        

#%%
class PAMSSite:
    def __init__(self, sitename):
        self.sitename = sitename
#%%
class FrontDetectorChromatogram(Chromatogram):
    pass

    
class BackDetectorChromatogram(Chromatogram):
    pass

#%%
class BaseSample:
    def __init__(self, front_path, back_path):
        self.front = FrontDetectorChromatogram(front_path)
        self.back = BackDetectorChromatogram(back_path)
# Specific sample types inherit from BaseSample
class Sample(BaseSample):
    pass

class Blank(BaseSample):
    pass

class CVS(BaseSample):
    pass

class RTS(BaseSample):
    pass

class LCS(BaseSample):
    pass
#%%
class Calibration(BaseSample):
    pass

class CalibrationCurve:
    def __init__(self):
        self.calruns= None
#%%
class DetectionLimit(BaseSample):
    pass

class DetectionLimitCombined:
    def __init__(self):
        self.detectionlimits = None
#%%
if __name__ == "__main__":
    c = pcp.Compound.from_cid(8058)
    print(c.molecular_formula)
    print(c.molecular_weight)
    print(c.elements)
    
    
    





