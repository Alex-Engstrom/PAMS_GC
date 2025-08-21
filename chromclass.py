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

class Chromatogram:
    """Handles chromatogram data from .cdf files."""

    def __init__(self, filename, dataformat="cdf"):
        self.format = dataformat
        self.filename = Path(filename)
        self._datetime = None
        self._chromatogram = None
        self._peakareas = None
        self._loaded = False
        self._error = None

    def load(self):
        """Explicitly load all data"""
        try:
            self._datetime = self._get_datetime()
            self._chromatogram = self._generate_chrom()
            self._peakareas = self._generate_areas()
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
    def peakareas(self):
        if not self._loaded and self._peakareas is None:
            self._peakareas = self._generate_areas()
        return self._peakareas

        
    def _get_datetime(self):
            """Private method to get datetime"""
            # Implement your datetime extraction logic here
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
                required_vars = ["ordinate_values", "actual_run_time_length", "actual_delay_time"]
                for var in required_vars:
                    if var not in rootgrp.variables:
                        raise KeyError(f"Missing '{var}' in {self.filename}")
                
                signal = np.array(rootgrp.variables["ordinate_values"][:])
                runtime = float(rootgrp.variables["actual_run_time_length"][0])
                starttime = float(rootgrp.variables["actual_delay_time"][:])
                
                rt = np.linspace(start=starttime, stop=runtime, num=len(signal))
                return np.vstack((rt, signal))
        except Exception as e:
            print(f"Error accessing file {self.filename}: {e}")
            return None

    def _generate_areas(self):
        """Private method to generate peak areas"""
        try:
            with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
                if "peak_name" not in rootgrp.variables or "peak_amount" not in rootgrp.variables:
                    return pd.DataFrame(columns=["peak_name", "peak_amount"])
                
                peak_names = [
                    b"".join(np.ma.filled(row, b"")).decode("utf-8").strip()
                    for row in rootgrp.variables["peak_name"][:]
                ]

                peak_amount = rootgrp.variables["peak_amount"][:]

                return pd.DataFrame({
                    "peak_name": peak_names,
                    "peak_amount": peak_amount
                })
        except Exception as e:
            print(f"Error reading peak areas: {e}")
            return pd.DataFrame(columns=["peak_name", "peak_amount"])
    def list_netcdf_variables(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return rootgrp.variables.keys()
    def list_netcdf_attributes(self):
        with ncdf.Dataset(self.filename, "r", format="NETCDF3_CLASSIC") as rootgrp:
            return rootgrp.ncattrs()

    def __repr__(self):
        return f"Chromatogram({self.filename.name})"
 
class Day:
    """Represents a single day of chromatogram data."""

    def __init__(self, site_path, date):
        """
        Creates Day object that contains chromatogram objects collected that day.
        Parameters:
            site_path: str | Path -> path to site root
            date: str -> date in the format YYYYMMDD
        """
        self.date = datetime.strptime(date, "%Y%m%d")
        self.folder = Path(site_path) / str(self.date.year) / f"{self.date.month:02d}" / f"{self.date.day:02d}"
        self.filename = self.folder.stem
        self.chromatograms = self.get_samples(self.folder)
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
     
    def get_samples(self, path):
        """
        Returns a list of Sample/Blank/RTS/CVS/LCS objects
        by pairing Front and Back Signal files.
        """
        if not path.exists():
            print(f"Warning: Folder does not exist → {path}")
            return []

        # Step 1: scan all files and parse metadata
        files_by_run = defaultdict(dict)
        for file in path.glob("*.cdf"):
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
            "m":Calibration
        }

        for run_key, paths in files_by_run.items():
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
    def __init__(self):
        self.sitename = None
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
    d1 = Day(r"C:\AutoGCData\RB", "20250802")
    print(f"Day folder: {d1}")
    print(f"Chromatograms found: {len(d1.chromatograms)}")
    for chrom in d1.chromatograms:
        if type(chrom) == Sample:
            print("SAMPLE", chrom.front.datetime, chrom.front.filename,chrom.back.datetime, chrom.back.filename)
        elif type(chrom) == CVS:
            print("CVS", chrom.front.filename, chrom.back.filename)
            print(chrom.front.list_netcdf_variables())
            print(chrom.front.list_netcdf_attributes())




