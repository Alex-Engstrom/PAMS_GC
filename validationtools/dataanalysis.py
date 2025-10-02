# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 14:13:47 2025

@author: aengstrom
"""
#%% imports
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
import pubchempy as pcp
from collections import defaultdict
#%% Reference Tables
mdls = {
    "Ethane": 0.078055655,
    "Ethylene": 0.127390085,
    "Propane": 0.10797866,
    "Propylene": 0.082203853,
    "Iso-butane": 0.097361678,
    "N-butane": 0.087842483,
    "Acetylene": 0.055647304,
    "Trans-2-butene": 0.043039742,
    "1-butene": 0.041851477,
    "Cis-2-butene": 0.029815072,
    "Cyclopentane": 0.04629037,
    "Iso-pentane": 0.064683983,
    "N-pentane": 0.039727116,
    "1,3-butadiene": 0.057630161,
    "Trans-2-pentene": 0.033336211,
    "1-pentene": 0.019727555,
    "Cis-2-pentene": 0.038408349,
    "2,2-dimethylbutane": 0.031532306,
    "2,3-dimethylbutane": 0.039956442,
    "2-methylpentane": 0.042617874,
    "3-methylpentane": 0.028312098,
    "Isoprene": 0.040727808,
    "2-methyl-1-pentene": 0.055647304,
    "1-hexene": 0.045751712,
    "N-hexane": 0.103401919,
    "Methylcyclopentane": 0.053189589,
    "2,4-dimethylpentane": 0.045059345,
    "Benzene": 0.043663749,
    "Cyclohexane": 0.192023364,
    "2-methylhexane": 0.082184757,
    "2,3-dimethylpentane": 0.087019574,
    "3-methylhexane": 0.064261312,
    "2,2,4-trimethylpentane": 0.055415033,
    "N-heptane": 0.096396244,
    "Methylcyclohexane": 0.047209407,
    "2,3,4-trimethylpentane": 0.098142236,
    "Toluene": 0.063267951,
    "2-methylheptane": 0.082245563,
    "3-methylheptane": 0.087970654,
    "N-octane": 0.06676112,
    "Ethylbenzene": 0.043888929,
    "M&p-xylene": 0.11422051,
    "Styrene": 0.071308108,
    "O-xylene": 0.03827329,
    "N-nonane": 0.047331881,
    "Iso-propylbenzene": 0.031818511,
    "Alpha-pinene": 0.290982689,
    "N-propylbenzene": 0.049961085,
    "M-ethyltoluene": 0.1098,
    "P-ethyltoluene": 0.06758161,
    "1,3,5-tri-m-benzene": 0.064181312,
    "O-ethyltoluene": 0.050875593,
    "Beta-pinene": 0.094498014,
    "1,2,4-tri-m-benzene": 0.090868478,
    "N-decane": 0.043351442,
    "1,2,3-tri-m-benzene": 0.0623,
    "M-diethylbenzene": 0.047420049,
    "P-diethylbenzene": 0.0663,
    "N-undecane": 0.120920325,
    "N-dodecane": 0.107645674,
    "TNMTC": 10,
    "TNMHC": 10
}
aqs_file_guide = {"RD": 
                      [
    "Transaction Type",
    "Action Indicator",
    "State Code / Tribal Indicator",
    "County Code / Tribal Code",
    "Site Number",
    "Parameter",
    "POC",
    "Duration Code",
    "Reported Unit",
    "Method Code",
    "Sample Date",
    "Sample Begin Time",
    "Reported Sample Value",
    "Null Data Code",
    "Collection Frequency Code",
    "Monitor Protocol ID",
    "Qualifier Code - 1",
    "Qualifier Code - 2",
    "Qualifier Code - 3",
    "Qualifier Code - 4",
    "Qualifier Code - 5",
    "Qualifier Code - 6",
    "Qualifier Code - 7",
    "Qualifier Code - 8",
    "Qualifier Code - 9",
    "Qualifier Code - 10",
    "Alternate Method Detection Limit",
    "Uncertainty Value"
                        ],
    "RB": 
                        [
    "Transaction Type",
    "Action Indicator",
    "State Code / Tribal Indicator",
    "County Code / Tribal Code",
    "Site Number",
    "Parameter",
    "POC",
    "Duration Code",
    "Reported Unit",
    "Method Code",
    "Blank Type",
    "Blank Date",
    "Blank Time",
    "Blank Value",
    "Null Data Code",
    "Qualifier Code - 1",
    "Qualifier Code - 2",
    "Qualifier Code - 3",
    "Qualifier Code - 4",
    "Qualifier Code - 5",
    "Qualifier Code - 6",
    "Qualifier Code - 7",
    "Qualifier Code - 8",
    "Qualifier Code - 9",
    "Qualifier Code - 10",
    "Alternate Method Detection Limit",
    "Uncertainty Value"
                        ]
                }
aqs_compound_codes = {
    "Ethane": 43202,
    "Ethylene": 43203,
    "Propane": 43204,
    "Propylene": 43205,
    "Iso-butane": 43214,
    "N-butane": 43212,
    "Acetylene": 43206,
    "Trans-2-butene": 43216,
    "1-butene": 43280,
    "Cis-2-butene": 43217,
    "Cyclopentane": 43242,
    "Iso-pentane": 43221,
    "N-pentane": 43220,
    "1,3-butadiene": 43218,
    "Trans-2-pentene": 43226,
    "1-pentene": 43224,
    "Cis-2-pentene": 43227,
    "2,2-dimethylbutane": 43244,
    "2,3-dimethylbutane": 43284,
    "2-methylpentane": 43285,
    "3-methylpentane": 43230,
    "Isoprene": 43243,
    "2-methyl-1-pentene": 43246,
    "1-hexene": 43245,
    "N-hexane": 43231,
    "Methylcyclopentane": 43262,
    "2,4-dimethylpentane": 43247,
    "Benzene": 45201,
    "Cyclohexane": 43248,
    "2-methylhexane": 43263,
    "2,3-dimethylpentane": 43291,
    "3-methylhexane": 43249,
    "2,2,4-trimethylpentane": 43250,
    "N-heptane": 43232,
    "Methylcyclohexane": 43261,
    "2,3,4-trimethylpentane": 43252,
    "Toluene": 45202,
    "2-methylheptane": 43960,
    "3-methylheptane": 43253,
    "N-octane": 43233,
    "Ethylbenzene": 45203,
    "M&p-xylene": 45109,
    "Styrene": 45220,
    "O-xylene": 45204,
    "N-nonane": 43235,
    "Iso-propylbenzene": 45210,
    "Alpha-pinene": 43256,
    "N-propylbenzene": 45209,
    "M-ethyltoluene": 45212,
    "P-ethyltoluene": 45213,
    "1,3,5-tri-m-benzene": 45207,
    "O-ethyltoluene": 45211,
    "Beta-pinene": 43257,
    "1,2,4-tri-m-benzene": 45208,
    "N-decane": 43238,
    "1,2,3-tri-m-benzene": 45225,
    "M-diethylbenzene": 45218,
    "P-diethylbenzene": 45219,
    "N-undecane": 43954,
    "N-dodecane": 43141,
    "TNMTC": 43102,
    "TNMHC": 43000
}
plot_compounds = [
    "Ethane",
    "Ethylene",
    "Propane",
    "Propylene",
    "Iso-butane",
    "n-Butane",
    "Acetylene",
    "trans-2-Butene",
    "1-Butene",
    "cis-2-Butene",
    "Cyclopentane",
    "Iso-pentane",
    "n-Pentane",
    "1,3-Butadiene",
    "trans-2-Pentene",
    "1-Pentene",
    "cis-2-Pentene",
    "2,2-Dimethylbutane",
    "2,3-Dimethylbutane",
    "2-Methylpentane",
    "3-Methylpentane",
    "Isoprene",
    #"2-Methyl-1-Pentene",    #Not reported to AQS
    "1-Hexene"
    ]
bp_compounds = [
        "n-Hexane",
        "Methylcyclopentane",
        "2,4-Dimethylpentane",
        "Benzene",
        "Cyclohexane",
        "2-Methylhexane",
        "2,3-Dimethylpentane",
        "3-Methylhexane",
        "2,2,4-Trimethylpentane",
        "n-Heptane",
        "Methylcyclohexane",
        "2,3,4-Trimethylpentane",
        "Toluene",
        "2-Methylheptane",
        "3-Methylheptane",
        "n-Octane",
        "Ethylbenzene",
        "m&p-Xylene",
        "Styrene",
        "o-Xylene",
        "n-Nonane",
        "Iso-propylbenzene",
        #"alpha-Pinene",       #Not reported to AQS
        "n-Propylbenzene",
        "m-ethyltoluene",
        "p-Ethyltoluene",
        "1,3,5-Tri-m-benzene",
        "o-Ethyltoluene",
        #"beta-Pinene",     #Not reported to AQS
        "1,2,4-Tri-m-benzene",
        "n-Decane",
        "1,2,3-Tri-m-benzene",
        "m-Diethylbenzene",
        "p-Diethylbenzene",
        "n-Undecane",
        "n-Dodecane"
    ]
def reverse_aqs_compound_codes(aqs_code_dict = {
    "Ethane": 43202,
    "Ethylene": 43203,
    "Propane": 43204,
    "Propylene": 43205,
    "Iso-butane": 43214,
    "N-butane": 43212,
    "Acetylene": 43206,
    "Trans-2-butene": 43216,
    "1-butene": 43280,
    "Cis-2-butene": 43217,
    "Cyclopentane": 43242,
    "Iso-pentane": 43221,
    "N-pentane": 43220,
    "1,3-butadiene": 43218,
    "Trans-2-pentene": 43226,
    "1-pentene": 43224,
    "Cis-2-pentene": 43227,
    "2,2-dimethylbutane": 43244,
    "2,3-dimethylbutane": 43284,
    "2-methylpentane": 43285,
    "3-methylpentane": 43230,
    "Isoprene": 43243,
    "2-methyl-1-pentene": 43246,
    "1-hexene": 43245,
    "N-hexane": 43231,
    "Methylcyclopentane": 43262,
    "2,4-dimethylpentane": 43247,
    "Benzene": 45201,
    "Cyclohexane": 43248,
    "2-methylhexane": 43263,
    "2,3-dimethylpentane": 43291,
    "3-methylhexane": 43249,
    "2,2,4-trimethylpentane": 43250,
    "N-heptane": 43232,
    "Methylcyclohexane": 43261,
    "2,3,4-trimethylpentane": 43252,
    "Toluene": 45202,
    "2-methylheptane": 43960,
    "3-methylheptane": 43253,
    "N-octane": 43233,
    "Ethylbenzene": 45203,
    "M&p-xylene": 45109,
    "Styrene": 45220,
    "O-xylene": 45204,
    "N-nonane": 43235,
    "Iso-propylbenzene": 45210,
    "Alpha-pinene": 43256,
    "N-propylbenzene": 45209,
    "M-ethyltoluene": 45212,
    "P-ethyltoluene": 45213,
    "1,3,5-tri-m-benzene": 45207,
    "O-ethyltoluene": 45211,
    "Beta-pinene": 43257,
    "1,2,4-tri-m-benzene": 45208,
    "N-decane": 43238,
    "1,2,3-tri-m-benzene": 45225,
    "M-diethylbenzene": 45218,
    "P-diethylbenzene": 45219,
    "N-undecane": 43954,
    "N-dodecane": 43141,
    "TNMTC": 43102,
    "TNMHC": 43000
    }):
    return {value : key for key, value in aqs_code_dict.items()}
def plot_compound_list_generator():
    plot_cid_list = [
    6324,  # Ethane
    6325,  # Ethylene
    6334,  # Propane
    8252,  # Propylene
    6360,  # Iso-butane
    7843,  # n-Butane
    6326,  # Acetylene
    62695,  # trans-2-Butene
    7844,  # 1-Butene
    5287573,  # cis-2-Butene
    9253,  # Cyclopentane
    6556,  # Iso-pentane
    8003,  # n-Pentane
    7845,  # 1,3-Butadiene
    5326161,  # trans-2-Pentene
    8004,  # 1-Pentene
    5326160,  # cis-2-Pentene
    6403,  # 2,2-Dimethylbutane
    6589, # 2,3-Dimethylbutane 
    7892,  # 2-Methylpentane
    7282,  # 3-Methylpentane
    6557,  # Isoprene
    12986, # 2-Methyl-1-pentene
    11597  # 1-Hexene 
]
 
    compounds = [
    "Ethane",
    "Ethylene",
    "Propane",
    "Propylene",
    "Iso-butane",
    "n-Butane",
    "Acetylene",
    "trans-2-Butene",
    "1-Butene",
    "cis-2-Butene",
    "Cyclopentane",
    "Iso-pentane",
    "n-Pentane",
    "1,3-Butadiene",
    "trans-2-Pentene",
    "1-Pentene",
    "cis-2-Pentene",
    "2,2-Dimethylbutane",
    "2,3-Dimethylbutane",
    "2-Methylpentane",
    "3-Methylpentane",
    "Isoprene",
    "2-Methyl-1-Pentene",
    "1-Hexene"
    ]
    voc_cid = zip(compounds, plot_cid_list)
    compound_list = []
    for voc, cid in voc_cid:
        c = pcp.Compound.from_cid(cid)
        compound_dict = {voc :
                      {'molecular_formula': c.molecular_formula,
                                       'molecular_weight': c.molecular_weight,
                                       'carbon_content': c.elements.count('C'),
                                       'hydrogen_content': c.elements.count('H'),
                                       'cid': cid}
                      }
        compound_list.append(compound_dict)
        print(c.iupac_name)
    return compound_list
def bp_compound_list_generator():
    bp_cid_list = [8058, 7296, 7907, 241, 8078, 11582, 11260, 11507, 10907, 
                   8900, 7962, 11269, 1140, 11594, 11519, 356, 7500, 7929, 7501, 7237, 8141, 7406, 6654, 
                   7668, 12100, 12160, 7947, 11903, 14896, 7247, 15600, 10686, 8864, 7734, 14257, 8182]
 
    compounds = [
        "n-Hexane",
        "Methylcyclopentane",
        "2,4-Dimethylpentane",
        "Benzene",
        "Cyclohexane",
        "2-Methylhexane",
        "2,3-Dimethylpentane",
        "3-Methylhexane",
        "2,2,4-Trimethylpentane",
        "n-Heptane",
        "Methylcyclohexane",
        "2,3,4-Trimethylpentane",
        "Toluene",
        "2-Methylheptane",
        "3-Methylheptane",
        "n-Octane",
        "Ethylbenzene",
        "m&p-Xylene",
        "Styrene",
        "o-Xylene",
        "n-Nonane",
        "Iso-propylbenzene",
        "alpha-Pinene",
        "n-Propylbenzene",
        "m-ethyltoluene",
        "p-Ethyltoluene",
        "1,3,5-Tri-m-benzene",
        "o-Ethyltoluene",
        "beta-Pinene",
        "1,2,4-Tri-m-benzene",
        "n-Decane",
        "1,2,3-Tri-m-benzene",
        "m-Diethylbenzene",
        "p-Diethylbenzene",
        "n-Undecane",
        "n-Dodecane"
    ]
    voc_cid = zip(compounds, bp_cid_list)
    compound_list = []
    for voc, cid in voc_cid:
        c = pcp.Compound.from_cid(cid)
        compound_dict = {voc :
                      {'molecular_formula': c.molecular_formula,
                                       'molecular_weight': c.molecular_weight,
                                       'carbon_content': c.elements.count('C'),
                                       'hydrogen_content': c.elements.count('H'),
                                       'cid': cid}
                      }
        compound_list.append(compound_dict)
    return compound_list
#%% Class Definitions
class Chromatogram:
    """Handles chromatogram data from .cdf files."""
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
    def generate_blank_summary(self, mdls):
        mdl_df = pd.DataFrame.from_dict(mdls, orient = 'columns')
        df_columns = ['date_time']
        df_columns.extend([str(list(voc.keys())[0]) for voc in FrontDetectorChromatogram.plot_vocs])
        df_columns.extend([str(list(voc.keys())[0]) for voc in BackDetectorChromatogram.bp_vocs])
        blank_df = pd.DataFrame(columns = df_columns)
        for chrom in self.blank:
            front, back = chrom.front, chrom.back
            date_front, date_back = front.datetime, back.datetime
            if date_front != date_back:
                print('front and back chromatograms do not have the same date values') 
            else:
                voc_amounts = pd.concat([front.peakamounts, back.peakamounts], ignore_index =True)
                new_row = {col: None for col in blank_df.columns}
                new_row['date_time'] = date_front
                for index, row in voc_amounts.iterrows():
                    compound_name = row['peak_name']
                    amount = row['peak_amount']
                    if compound_name in blank_df.columns:
                        new_row[compound_name]= amount
                blank_df = pd.concat([blank_df, pd.DataFrame([new_row])], ignore_index = True)
        return blank_df
                    
                    

                        

                
                
                
                
                
            
        
            
        
    def __str__(self):
        return f"{self.folder}"
        

    
class PAMSSite(DataAnalysis):
    def __init__(self):
        self.test = None

class RB(DataAnalysis):     
    sitename = 'RB'
    aqscode = None       
    def __init__(self, site_path, year, month = None, day = None):
        self.mdl = self.get_mdl()
    def get_mdl(self):
        return pd.read_csv('rb_mdl.csv')

class FrontDetectorChromatogram(Chromatogram):
    rf = 6104
    plot_vocs = [
                  {'Ethane': {'molecular_formula': 'C2H6', 'molecular_weight': 30.07, 'carbon_content': 2, 'hydrogen_content': 6, 'cid': 6324}}, 
                  {'Ethylene': {'molecular_formula': 'C2H4', 'molecular_weight': 28.05, 'carbon_content': 2, 'hydrogen_content': 4, 'cid': 6325}}, 
                  {'Propane': {'molecular_formula': 'C3H8', 'molecular_weight': 44.1, 'carbon_content': 3, 'hydrogen_content': 8, 'cid': 6334}}, 
                  {'Propylene': {'molecular_formula': 'C3H6', 'molecular_weight': 42.08, 'carbon_content': 3, 'hydrogen_content': 6, 'cid': 8252}}, 
                  {'Iso-butane': {'molecular_formula': 'C4H10', 'molecular_weight': 58.12, 'carbon_content': 4, 'hydrogen_content': 10, 'cid': 6360}}, 
                  {'n-Butane': {'molecular_formula': 'C4H10', 'molecular_weight': 58.12, 'carbon_content': 4, 'hydrogen_content': 10, 'cid': 7843}}, 
                  {'Acetylene': {'molecular_formula': 'C2H2', 'molecular_weight': 26.04, 'carbon_content': 2, 'hydrogen_content': 2, 'cid': 6326}}, 
                  {'trans-2-Butene': {'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 62695}}, 
                  {'1-Butene': {'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 7844}}, 
                  {'cis-2-Butene': {'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 5287573}}, 
                  {'Cyclopentane': {'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 9253}}, 
                  {'Iso-pentane': {'molecular_formula': 'C5H12', 'molecular_weight': 72.15, 'carbon_content': 5, 'hydrogen_content': 12, 'cid': 6556}},
                  {'n-Pentane': {'molecular_formula': 'C5H12', 'molecular_weight': 72.15, 'carbon_content': 5, 'hydrogen_content': 12, 'cid': 8003}}, 
                  {'1,3-Butadiene': {'molecular_formula': 'C4H6', 'molecular_weight': 54.09, 'carbon_content': 4, 'hydrogen_content': 6, 'cid': 7845}}, 
                  {'trans-2-Pentene': {'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 5326161}},
                  {'1-Pentene': {'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 8004}}, 
                  {'cis-2-Pentene': {'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 5326160}}, 
                  {'2,2-Dimethylbutane': {'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 6403}}, 
                  {'2,3-Dimethylbutane': {'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 6589}}, 
                  {'2-Methylpentane': {'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 7892}}, 
                  {'3-Methylpentane': {'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 7282}},
                  {'Isoprene': {'molecular_formula': 'C5H8', 'molecular_weight': 68.12, 'carbon_content': 5, 'hydrogen_content': 8, 'cid': 6557}},
                  {'2-Methyl-1-Pentene': {'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 12986}}, 
                  {'1-Hexene': {'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 11597}}
                  ]
    pass

    
class BackDetectorChromatogram(Chromatogram):
    rf = 5828
    bp_vocs = [
               {'n-Hexane': {'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 8058}}, 
               {'Methylcyclopentane': {'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 7296}},
               {'2,4-Dimethylpentane': {'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 7907}}, 
               {'Benzene': {'molecular_formula': 'C6H6', 'molecular_weight': 78.11, 'carbon_content': 6, 'hydrogen_content': 6, 'cid': 241}}, 
               {'Cyclohexane': {'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 8078}}, 
               {'2-Methylhexane': {'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11582}}, 
               {'2,3-Dimethylpentane': {'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11260}}, 
               {'3-Methylhexane': {'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11507}}, 
               {'2,2,4-Trimethylpentane': {'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 10907}}, 
               {'n-Heptane': {'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 8900}}, 
               {'Methylcyclohexane': {'molecular_formula': 'C7H14', 'molecular_weight': 98.19, 'carbon_content': 7, 'hydrogen_content': 14, 'cid': 7962}}, 
               {'2,3,4-Trimethylpentane': {'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11269}},
               {'Toluene': {'molecular_formula': 'C7H8', 'molecular_weight': 92.14, 'carbon_content': 7, 'hydrogen_content': 8, 'cid': 1140}}, 
               {'2-Methylheptane': {'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11594}},
               {'3-Methylheptane': {'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11519}}, 
               {'n-Octane': {'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 356}}, 
               {'Ethylbenzene': {'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7500}}, 
               {'m&p-Xylene': {'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7929}}, 
               {'Styrene': {'molecular_formula': 'C8H8', 'molecular_weight': 104.15, 'carbon_content': 8, 'hydrogen_content': 8, 'cid': 7501}},
               {'o-Xylene': {'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7237}}, 
               {'n-Nonane': {'molecular_formula': 'C9H20', 'molecular_weight': 128.25, 'carbon_content': 9, 'hydrogen_content': 20, 'cid': 8141}}, 
               {'Iso-propylbenzene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7406}}, 
               {'alpha-Pinene': {'molecular_formula': 'C10H16', 'molecular_weight': 136.23, 'carbon_content': 10, 'hydrogen_content': 16, 'cid': 6654}},
               {'n-Propylbenzene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7668}}, 
               {'m-ethyltoluene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 12100}},
               {'p-Ethyltoluene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 12160}}, 
               {'1,3,5-Tri-m-benzene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7947}}, 
               {'o-Ethyltoluene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 11903}}, 
               {'beta-Pinene': {'molecular_formula': 'C10H16', 'molecular_weight': 136.23, 'carbon_content': 10, 'hydrogen_content': 16, 'cid': 14896}}, 
               {'1,2,4-Tri-m-benzene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7247}}, 
               {'n-Decane': {'molecular_formula': 'C10H22', 'molecular_weight': 142.28, 'carbon_content': 10, 'hydrogen_content': 22, 'cid': 15600}}, 
               {'1,2,3-Tri-m-benzene': {'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 10686}},
               {'m-Diethylbenzene': {'molecular_formula': 'C10H14', 'molecular_weight': 134.22, 'carbon_content': 10, 'hydrogen_content': 14, 'cid': 8864}}, 
               {'p-Diethylbenzene': {'molecular_formula': 'C10H14', 'molecular_weight': 134.22, 'carbon_content': 10, 'hydrogen_content': 14, 'cid': 7734}},
               {'n-Undecane': {'molecular_formula': 'C11H24', 'molecular_weight': 156.31, 'carbon_content': 11, 'hydrogen_content': 24, 'cid': 14257}}, 
               {'n-Dodecane': {'molecular_formula': 'C12H26', 'molecular_weight': 170.33, 'carbon_content': 12, 'hydrogen_content': 26, 'cid': 8182}}
               ]
    pass


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

class Calibration(BaseSample):
    pass

class CalibrationCurve:
    def __init__(self):
        self.calruns= None
        
class DetectionLimit(BaseSample):
    pass

class DetectionLimitCombined:
    def __init__(self):
        self.detectionlimits = None

#%% Data Analysis Functions
def generate_df_from_aqs_file(aqs_file, transaction_type: str):
    column_names = aqs_file_guide[transaction_type]
    return pd.read_csv(aqs_file, names = column_names, index_col = False, sep = "|", dtype = 'str')
def summarize_blanks(aqs_file):
    aqs_df = generate_df_from_aqs_file(aqs_file=aqs_file, transaction_type="RB")
    aqs_df.columns = aqs_df.columns.str.replace(' ', '_').str.lower()
    aqs_df['blank_date'] = pd.to_datetime(aqs_df['blank_date'], format='%Y%m%d')
    aqs_df['blank_date'] = aqs_df['blank_date'].apply(str)
    qualifier_columns = [col for col in aqs_df.columns if col.startswith('qualifier_code')]
    
    # Filter rows where any qualifier column contains 'LB'
    lb_mask = (aqs_df[qualifier_columns] == 'LB').any(axis=1)
    lb_df = aqs_df[lb_mask]
    aqs_code_to_voc_name = reverse_aqs_compound_codes()
    # Group by blank_date and collect parameters
    date_series = lb_df.groupby('blank_date')['parameter'].apply(lambda x: [aqs_code_to_voc_name.get(int(p)) for p in x]).reset_index(name = 'failing_blanks')
    # Returns a series with an the index as the datetime object and the value a list of vocs with LB flags on that day
    return date_series

def summarize_nulls(aqs_file):
    aqs_df = generate_df_from_aqs_file(aqs_file=aqs_file, transaction_type="RD")
    aqs_df.columns = aqs_df.columns.str.replace(' ', '_').str.lower()
    aqs_df['sample_datetime'] = aqs_df['sample_date']+' '+aqs_df['sample_begin_time']
    aqs_df['sample_datetime'] = pd.to_datetime(aqs_df['sample_datetime'], errors='coerce')
    needed_columns = ['sample_datetime', 'null_data_code', 'parameter']
    aqs_df_clean = aqs_df[needed_columns]
    by_null_code = aqs_df_clean.groupby(['null_data_code', 'sample_datetime']).agg(list).reset_index()
    nulls_by_date = aqs_df_clean.groupby(['sample_datetime'])['null_data_code'].count().reset_index()
    nulls_by_count = nulls_by_date.groupby(['null_data_code'])['sample_datetime'].agg(list).reset_index()
    full_hour_nulls = {}
    column_nulls = {'plot': {}, 'bp': {}}
    compound_nulls = {}
    for code in [code for code in by_null_code['null_data_code'].unique() if code not in ['AY', 'TC']]:
        mask = by_null_code['null_data_code'] == code
        masked_by_null_code = by_null_code[mask]
        for index,row in masked_by_null_code.iterrows():
            date_time = row['sample_datetime']
            compounds = row['parameter']
            
            if len(compounds) == 58:
                if code not in full_hour_nulls:
                    full_hour_nulls[code] = []
                full_hour_nulls[code].append(date_time)
                
            elif len(compounds) < 58:
                aqs_code_to_voc_name = reverse_aqs_compound_codes()
                compound_names = [aqs_code_to_voc_name.get(int(p)) for p in compounds if p not in ['43000','43102']]
                # Make sure all names are strings and handle None values
                compound_names = [name for name in compound_names if name is not None]
                
                if sorted([name.lower() for name in compound_names]) == sorted([name.lower() for name in plot_compounds]):
                    if code not in column_nulls['plot']:
                        column_nulls['plot'][code] = []
                    column_nulls['plot'][code].append(date_time)
                elif sorted([name.lower() for name in compound_names]) == sorted([name.lower() for name in bp_compounds]):
                    if code not in column_nulls['bp']:
                        column_nulls['bp'][code] = []
                    column_nulls['bp'][code].append(date_time)
                else:
                    compound_names = [aqs_code_to_voc_name.get(int(p)) for p in compounds]
                    print(date_time, compound_names)
                
    return  full_hour_nulls, nulls_by_count

def summarize_qualifiers(aqs_file):
    aqs_df = generate_df_from_aqs_file(aqs_file=aqs_file, transaction_type="RD")
    aqs_df.columns = aqs_df.columns.str.replace(' ', '_').str.lower()
    aqs_df['sample_datetime'] = aqs_df['sample_date']+' '+aqs_df['sample_begin_time']
    aqs_df['sample_datetime'] = pd.to_datetime(aqs_df['sample_datetime'], errors='coerce')
    
    # Melt the code columns into a long format
    code_cols = [col for col in aqs_df.columns if col.startswith("qualifier")]
    df_long = aqs_df.melt(
        id_vars=['sample_datetime', 'parameter'],
        value_vars=code_cols,
        value_name='code'
    ).dropna(subset=['code'])
    
    # Now we can group by compound and code
    summary = (
        df_long
        .groupby(['parameter', 'code'])
        .agg(
            count=('code', 'count'),          # How many times this code was applied
            dates=('sample_datetime', lambda x: list(x)) # List of dates
        )
        .reset_index()
    )
    mask = ~summary['code'].isin(['MD', 'SQ', 'ND', 'LB'])
    return summary[mask]
    
# #%% Main
# if __name__ == "__main__":
#     #full_hour_nulls, nulls_by_count = summarize_nulls(r"C:\Users\aengstrom\Desktop\RD_SITE-RB_RedButte_INSERT_08-01-2025_to_08-31-2025.txt")
#     #blank_dict = summarize_blanks(aqs_file=r"C:\Users\aengstrom\Desktop\RB_SITE-RB_RedButte_INSERT_08-01-2025_to_08-31-2025.txt")
#     aqs_file_rd = r"C:\Users\aengstrom\Desktop\RD_SITE-RB_RedButte_INSERT_08-01-2025_to_08-31-2025.txt"
#     code_to_name = dt.reverse_aqs_compound_codes()
#     qualifiers = dt.summarize_qualifiers(aqs_file_rd)
#     for code in qualifiers['code'].unique():
#         print(f"- **{code}**")
#         mask = qualifiers['code'] == code
#         masked = qualifiers[mask]
#         for index, row in masked.iterrows():
#             compound_code = int(row['parameter'])
#             hours = row['count']
#             voc_name =  code_to_name.get(compound_code)
#             print(f"  - {voc_name} : {hours} hour(s) qualified")



    

        

    
    
    





