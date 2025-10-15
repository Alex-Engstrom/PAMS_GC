
import netCDF4 as ncdf
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
import re
from collections import defaultdict
from dateutil import parser

#%% Class Definitions
class Chromatogram:
    aqs_compound_to_code = {
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
        "TNMHC": 43000,
        'Bp unid': 10000,
        'Plot unid':  20000
    }
    aqs_code_to_compound = {
                                   43202: 'Ethane', 
                                   43203: 'Ethylene', 
                                   43204: 'Propane', 
                                   43205: 'Propylene', 
                                   43214: 'Iso-butane', 
                                   43212: 'N-butane', 
                                   43206: 'Acetylene', 
                                   43216: 'Trans-2-butene', 
                                   43280: '1-butene', 
                                   43217: 'Cis-2-butene', 
                                   43242: 'Cyclopentane', 
                                   43221: 'Iso-pentane', 
                                   43220: 'N-pentane', 
                                   43218: '1,3-butadiene', 
                                   43226: 'Trans-2-pentene', 
                                   43224: '1-pentene', 
                                   43227: 'Cis-2-pentene', 
                                   43244: '2,2-dimethylbutane', 
                                   43284: '2,3-dimethylbutane',
                                   43285: '2-methylpentane', 
                                   43230: '3-methylpentane', 
                                   43243: 'Isoprene',
                                   43246: '2-methyl-1-pentene', 
                                   43245: '1-hexene', 
                                   43231: 'N-hexane', 
                                   43262: 'Methylcyclopentane', 
                                   43247: '2,4-dimethylpentane', 
                                   45201: 'Benzene',
                                   43248: 'Cyclohexane', 
                                   43263: '2-methylhexane',
                                   43291: '2,3-dimethylpentane', 
                                   43249: '3-methylhexane', 
                                   43250: '2,2,4-trimethylpentane', 
                                   43232: 'N-heptane', 
                                   43261: 'Methylcyclohexane', 
                                   43252: '2,3,4-trimethylpentane',
                                   45202: 'Toluene', 
                                   43960: '2-methylheptane', 
                                   43253: '3-methylheptane', 
                                   43233: 'N-octane', 
                                   45203: 'Ethylbenzene', 
                                   45109: 'M&p-xylene', 
                                   45220: 'Styrene', 
                                   45204: 'O-xylene', 
                                   43235: 'N-nonane', 
                                   45210: 'Iso-propylbenzene', 
                                   43256: 'Alpha-pinene', 
                                   45209: 'N-propylbenzene', 
                                   45212: 'M-ethyltoluene', 
                                   45213: 'P-ethyltoluene', 
                                   45207: '1,3,5-tri-m-benzene',
                                   45211: 'O-ethyltoluene', 
                                   43257: 'Beta-pinene', 
                                   45208: '1,2,4-tri-m-benzene', 
                                   43238: 'N-decane', 
                                   45225: '1,2,3-tri-m-benzene', 
                                   45218: 'M-diethylbenzene',
                                   45219: 'P-diethylbenzene', 
                                   43954: 'N-undecane', 
                                   43141: 'N-dodecane', 
                                   43102: 'TNMTC', 
                                   43000: 'TNMHC',
                                   10000: 'BP UnID',
                                   20000: 'PLOT UnID'}
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
                df = pd.DataFrame(data)
                mask = df['peak_name'].str.len() > 0
                df = df[mask]
                df['peak_name'] = df['peak_name'].map(lambda x: Chromatogram.aqs_compound_to_code.get(x.capitalize()))
                if attribute == 'peakamounts':
                    mask = ~df['peak_name'].isin([10000,20000])
                    tnmtc = df[mask]['peak_amount'].sum()
                    tnmhc = df['peak_amount'].sum()
                    df.loc[len(df)] = [43102, tnmtc]
                    df.loc[len(df)] = [43000, tnmhc]
                    df = df.astype({'peak_name' : int})
                return df
                

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
    compound_categories = {
        "Ethane": "ALKANE",
        "Ethylene": "ALKENE",
        "Propane": "ALKANE",
        "Propylene": "ALKENE",
        "Iso-butane": "ALKANE",
        "N-butane": "ALKANE",
        "Acetylene": "ALKYNE",
        "Trans-2-butene": "ALKENE",
        "1-butene": "ALKENE",
        "Cis-2-butene": "ALKENE",
        "Cyclopentane": "ALKANE",
        "Iso-pentane": "ALKANE",
        "N-pentane": "ALKANE",
        "1,3-butadiene": "ALKENE",
        "Trans-2-pentene": "ALKENE",
        "1-pentene": "ALKENE",
        "Cis-2-pentene": "ALKENE",
        "2,2-dimethylbutane": "ALKANE",
        "2,3-dimethylbutane": "ALKANE",
        "2-methylpentane": "ALKANE",
        "3-methylpentane": "ALKANE",
        "Isoprene": "TERPENE",
        "2-methyl-1-pentene": "ALKENE",
        "1-hexene": "ALKENE",
        "N-hexane": "ALKANE",
        "Methylcyclopentane": "ALKANE",
        "2,4-dimethylpentane": "ALKANE",
        "Benzene": "AROMATIC",
        "Cyclohexane": "ALKANE",
        "2-methylhexane": "ALKANE",
        "2,3-dimethylpentane": "ALKANE",
        "3-methylhexane": "ALKANE",
        "2,2,4-trimethylpentane": "ALKANE",
        "N-heptane": "ALKANE",
        "Methylcyclohexane": "ALKANE",
        "2,3,4-trimethylpentane": "ALKANE",
        "Toluene": "AROMATIC",
        "2-methylheptane": "ALKANE",
        "3-methylheptane": "ALKANE",
        "N-octane": "ALKANE",
        "Ethylbenzene": "AROMATIC",
        "M&p-xylene": "AROMATIC",
        "Styrene": "AROMATIC",
        "O-xylene": "AROMATIC",
        "N-nonane": "ALKANE",
        "Iso-propylbenzene": "AROMATIC",
        "Alpha-pinene": "TERPENE",
        "N-propylbenzene": "AROMATIC",
        "M-ethyltoluene": "AROMATIC",
        "P-ethyltoluene": "AROMATIC",
        "1,3,5-tri-m-benzene": "AROMATIC",
        "O-ethyltoluene": "AROMATIC",
        "Beta-pinene": "TERPENE",
        "1,2,4-tri-m-benzene": "AROMATIC",
        "N-decane": "ALKANE",
        "1,2,3-tri-m-benzene": "AROMATIC",
        "M-diethylbenzene": "AROMATIC",
        "P-diethylbenzene": "AROMATIC",
        "N-undecane": "ALKANE",
        "N-dodecane": "ALKANE"
    }
    def __init__(self, site_name, data_path, year, month = None, day = None, use_filesystem = False, use_compound_code = False):
        """
        Creates DataAnalysis object that contains chromatogram objects collected over a year, month, or day.
        Parameters:
            data_path: str | Path -> path to site root
            year: str or int-> year of interest. If no other date inputs are entered, data from the whole year will be extracted
            month: str or int -> month of interest.
            day: str or int -> day of interest
        """
        self.use_compound_code = use_compound_code
        self.use_filesystem = use_filesystem
        self.date = [year, month, day]
        self.folder = self._find_root_dir(data_path = data_path, date = self.date)
        self.filename = self.folder.stem
        self.ambient = self._get_chromatograms(self.folder, 's')
        self.cvs = self._get_chromatograms(self.folder, 'c')
        self.lcs = self._get_chromatograms(self.folder, 'e')
        self.blank = self._get_chromatograms(self.folder, 'b')
        self.rts = self._get_chromatograms(self.folder, 'q')
        
    def _find_root_dir(self, data_path, date: list):
        if self.use_filesystem:
            root_dir = Path(data_path) / str(date[0])
            for j in date[1:]:
                if j:
                    root_dir = root_dir / f"{int(j):02d}"
        else:
            root_dir = Path(data_path)
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
            "s": Ambient,
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

    def generate_data_summary(self, chrom_type, mdls):
        mdl_df = pd.DataFrame(mdls, index=[0])
        mdl_df.columns = mdl_df.columns.map(str)
        mdl_df = mdl_df.rename(columns = lambda x: Chromatogram.aqs_compound_to_code.get(x,x))
        # Build column list
        df_columns = ['date_time']
        df_columns.extend(list(FrontDetectorChromatogram.plot_vocs.keys()))
        df_columns.extend(list(BackDetectorChromatogram.bp_vocs.keys()))
        chrom_type_guide = {'blank': self.blank,
                            'ambient': self.ambient,
                            'cvs': self.cvs,
                            'rts': self.rts,
                            'lcs': self.lcs
                            }
        # Process all chromatograms
        rows = []
        for chrom in chrom_type_guide[chrom_type]:
            if chrom.front.datetime != chrom.back.datetime:
                print('front and back chromatograms do not have the same date values')
                continue
                
            # Combine peak amounts and convert to dict for easier lookup
            voc_amounts = pd.concat([chrom.front.peakamounts, chrom.back.peakamounts], ignore_index=True)
            amount_dict = voc_amounts.set_index('peak_name')['peak_amount'].to_dict()
            
            # Build row
            row = {'date_time': chrom.front.datetime}
            row.update({col: amount_dict.get(col) for col in df_columns if col != 'date_time'})
            rows.append(row)
        
        blank_df = pd.DataFrame(rows, columns=df_columns)
        
        return blank_df, mdl_df
    
    def compounds_above_mdl(self, mdls):
        blank_df, mdl_df = self.generate_data_summary(chrom_type = 'blank', mdls = mdls)    
        compound_columns = [col for col in blank_df.columns if col != 'date_time']
        mdl_series = mdl_df.iloc[0]
        
        results = []
        for _, row in blank_df.iterrows():
            compounds_above = [
                compound for compound in compound_columns
                if (compound in mdl_series.index and 
                    pd.notna(row[compound]) and 
                    row[compound] > mdl_series[compound])
            ]
            results.append({
                'date_time': row['date_time'],
                'compounds_above_mdl': compounds_above,
                'count_above_mdl': len(compounds_above)
            })
        
        return pd.DataFrame(results)
    
    def sort_by_type(self):
        compound_categories = DataAnalysis.compound_categories
        compound_df = pd.DataFrame(list(compound_categories.items()), columns=['compound', 'category'])
        group_by_category = compound_df.groupby(['category'])['compound'].apply(lambda x: [Chromatogram.aqs_compound_to_code[c] for c in x])
        return group_by_category
    
    def check_ratios(self, mdls):
        voc_df, mdl_df = self.generate_data_summary(chrom_type = 'ambient', mdls = mdls)                                          
        voc_df.iloc[:,1:] = voc_df.iloc[:,1:].apply(pd.to_numeric)
        mdl_number = mdl_df.squeeze()       # convert 1-row DataFrame → Series
        mdl_number.index = mdl_number.index.map(int)  # ensure numeric index (AQS codes)
        compound_by_type = self.sort_by_type()
        alkanes = [compound for compound in compound_by_type['ALKANE']]
        alkenes = [compound for compound in compound_by_type['ALKENE']]
        # Commented out screens that do not require qualififcation per TAD table 10-1
        conditions = [
            ((voc_df[45201] > voc_df[45202]) & (voc_df[45201] > 3*mdl_number[45201]), 'benzene_gt_toluene',[45201,45202]),
            ((voc_df[45201] > voc_df[43202]) & (voc_df[45201] > 3*mdl_number[45201]), 'benzene_gt_ethane', [45201,43202]),
            ((voc_df[43203] > voc_df[43202]) & (voc_df[43203] > 3*mdl_number[43203]), 'ethylene_gt_ethane', [43203,43202]),
            ((voc_df[43205] > voc_df[43204]) & (voc_df[43205] > 3*mdl_number[43205]), 'propylene_gt_propane', [43205,43204]),
            ((voc_df[45204] > voc_df[45109]) & (voc_df[45204] > 3*mdl_number[45204]), 'oxylene_gt_mpxylene', [45204,45109]),
            #((voc_df[43263] < voc_df[43291]) & (voc_df[43291] > 3*mdl_number[43291]), '23dimethylpentane_gt_2methylhexane', [43263, 43291]),
            #((voc_df[43262] < voc_df[43247]) & (voc_df[43247] > 3*mdl_number[43247]), '24dimethylpentane_gt_methylcyclopentane', [43262, 43247]),
            #((voc_df[43214] > voc_df[43212]) & (voc_df[43214] > 3*mdl_number[43214]), 'isobutane_gt_nbutane', [43214, 43212]),
            ((voc_df[43230] > .6*voc_df[43285]) & (voc_df[43230] > 3*mdl_number[43230]), '3methylpentane_gt_2methylpentane', [43230,43285]),
            ((voc_df[43954] > voc_df[43238]) & (voc_df[43954] > 3*mdl_number[43954]), 'nundecane_gt_ndecane', [43954, 43238]),
            #(~((voc_df[43221] > voc_df[43220]) & (voc_df[43220] > voc_df[43242])) & (voc_df[43221] > 3*mdl_number[43221]) & (voc_df[43220] > 3*mdl_number[43220]) & (voc_df[43242] > 3*mdl_number[43242]), 
            # 'not_isopentane_gt_npentane_gt_cyclopentane', [43221,43220,43242]),
            (voc_df[alkenes].sum(axis = 1) > voc_df[alkanes].sum(axis = 1), 'alkenes_gt_alkanes', ['alkanes', 'alkenes']) 
             ]
        flagged = []
        for cond, label, compounds in conditions:
            subset = voc_df.loc[cond].copy()
            subset['screen_reason'] = label
            subset['compounds'] = [[Chromatogram.aqs_code_to_compound[compound] for compound in compounds if compound in Chromatogram.aqs_code_to_compound.keys()]]*len(subset)
            flagged.append(subset)

        ratios_check = pd.concat(flagged, ignore_index=True)
        ratios = ratios_check[['date_time', 'screen_reason','compounds']].copy()
        return ratios
  
    def __str__(self):
        return f"{self.folder}"
        

    
class PAMSSite(DataAnalysis):
    def __init__(self):
        self.test = None

class RB(DataAnalysis):     
    sitename = 'RB'
    aqscode = None       
    def __init__(self, data_path, year, month = None, day = None):
        self.mdl = self.get_mdl()
    def get_mdl(self):
        return pd.read_csv('rb_mdl.csv')

class FrontDetectorChromatogram(Chromatogram):
    
    plot_vocs = {
        43202: {'voc_name': 'Ethane', 'molecular_formula': 'C2H6', 'molecular_weight': 30.07, 'carbon_content': 2, 'hydrogen_content': 6, 'cid': 6324},
        43203: {'voc_name': 'Ethylene', 'molecular_formula': 'C2H4', 'molecular_weight': 28.05, 'carbon_content': 2, 'hydrogen_content': 4, 'cid': 6325},
        43204: {'voc_name': 'Propane', 'molecular_formula': 'C3H8', 'molecular_weight': 44.1, 'carbon_content': 3, 'hydrogen_content': 8, 'cid': 6334},
        43205: {'voc_name': 'Propylene', 'molecular_formula': 'C3H6', 'molecular_weight': 42.08, 'carbon_content': 3, 'hydrogen_content': 6, 'cid': 8252},
        43214: {'voc_name': 'Iso-butane', 'molecular_formula': 'C4H10', 'molecular_weight': 58.12, 'carbon_content': 4, 'hydrogen_content': 10, 'cid': 6360},
        43212: {'voc_name': 'N-butane', 'molecular_formula': 'C4H10', 'molecular_weight': 58.12, 'carbon_content': 4, 'hydrogen_content': 10, 'cid': 7843},
        43206: {'voc_name': 'Acetylene', 'molecular_formula': 'C2H2', 'molecular_weight': 26.04, 'carbon_content': 2, 'hydrogen_content': 2, 'cid': 6326},
        43216: {'voc_name': 'Trans-2-butene', 'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 62695},
        43280: {'voc_name': '1-butene', 'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 7844},
        43217: {'voc_name': 'Cis-2-butene', 'molecular_formula': 'C4H8', 'molecular_weight': 56.11, 'carbon_content': 4, 'hydrogen_content': 8, 'cid': 5287573},
        43242: {'voc_name': 'Cyclopentane', 'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 9253},
        43221: {'voc_name': 'Iso-pentane', 'molecular_formula': 'C5H12', 'molecular_weight': 72.15, 'carbon_content': 5, 'hydrogen_content': 12, 'cid': 6556},
        43220: {'voc_name': 'N-pentane', 'molecular_formula': 'C5H12', 'molecular_weight': 72.15, 'carbon_content': 5, 'hydrogen_content': 12, 'cid': 8003},
        43218: {'voc_name': '1,3-butadiene', 'molecular_formula': 'C4H6', 'molecular_weight': 54.09, 'carbon_content': 4, 'hydrogen_content': 6, 'cid': 7845},
        43226: {'voc_name': 'Trans-2-pentene', 'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 5326161},
        43224: {'voc_name': '1-pentene', 'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 8004},
        43227: {'voc_name': 'Cis-2-pentene', 'molecular_formula': 'C5H10', 'molecular_weight': 70.13, 'carbon_content': 5, 'hydrogen_content': 10, 'cid': 5326160},
        43244: {'voc_name': '2,2-dimethylbutane', 'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 6403},
        43284: {'voc_name': '2,3-dimethylbutane', 'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 6589},
        43285: {'voc_name': '2-methylpentane', 'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 7892},
        43230: {'voc_name': '3-methylpentane', 'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 7282},
        43243: {'voc_name': 'Isoprene', 'molecular_formula': 'C5H8', 'molecular_weight': 68.12, 'carbon_content': 5, 'hydrogen_content': 8, 'cid': 6557},
        43246: {'voc_name': '2-methyl-1-pentene', 'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 12986},
        43245: {'voc_name': '1-hexene', 'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 11597}
    }
    pass

    
class BackDetectorChromatogram(Chromatogram):

    bp_vocs = {
        43231: {'voc_name': 'N-hexane', 'molecular_formula': 'C6H14', 'molecular_weight': 86.18, 'carbon_content': 6, 'hydrogen_content': 14, 'cid': 8058},
        43262: {'voc_name': 'Methylcyclopentane', 'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 7296},
        43247: {'voc_name': '2,4-dimethylpentane', 'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 7907},
        45201: {'voc_name': 'Benzene', 'molecular_formula': 'C6H6', 'molecular_weight': 78.11, 'carbon_content': 6, 'hydrogen_content': 6, 'cid': 241},
        43248: {'voc_name': 'Cyclohexane', 'molecular_formula': 'C6H12', 'molecular_weight': 84.16, 'carbon_content': 6, 'hydrogen_content': 12, 'cid': 8078},
        43263: {'voc_name': '2-methylhexane', 'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11582},
        43291: {'voc_name': '2,3-dimethylpentane', 'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11260},
        43249: {'voc_name': '3-methylhexane', 'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 11507},
        43250: {'voc_name': '2,2,4-trimethylpentane', 'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 10907},
        43232: {'voc_name': 'N-heptane', 'molecular_formula': 'C7H16', 'molecular_weight': 100.2, 'carbon_content': 7, 'hydrogen_content': 16, 'cid': 8900},
        43261: {'voc_name': 'Methylcyclohexane', 'molecular_formula': 'C7H14', 'molecular_weight': 98.19, 'carbon_content': 7, 'hydrogen_content': 14, 'cid': 7962},
        43252: {'voc_name': '2,3,4-trimethylpentane', 'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11269},
        45202: {'voc_name': 'Toluene', 'molecular_formula': 'C7H8', 'molecular_weight': 92.14, 'carbon_content': 7, 'hydrogen_content': 8, 'cid': 1140},
        43960: {'voc_name': '2-methylheptane', 'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11594},
        43253: {'voc_name': '3-methylheptane', 'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 11519},
        43233: {'voc_name': 'N-octane', 'molecular_formula': 'C8H18', 'molecular_weight': 114.23, 'carbon_content': 8, 'hydrogen_content': 18, 'cid': 356},
        45203: {'voc_name': 'Ethylbenzene', 'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7500},
        45109: {'voc_name': 'M&p-xylene', 'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7929},
        45220: {'voc_name': 'Styrene', 'molecular_formula': 'C8H8', 'molecular_weight': 104.15, 'carbon_content': 8, 'hydrogen_content': 8, 'cid': 7501},
        45204: {'voc_name': 'O-xylene', 'molecular_formula': 'C8H10', 'molecular_weight': 106.16, 'carbon_content': 8, 'hydrogen_content': 10, 'cid': 7237},
        43235: {'voc_name': 'N-nonane', 'molecular_formula': 'C9H20', 'molecular_weight': 128.25, 'carbon_content': 9, 'hydrogen_content': 20, 'cid': 8141},
        45210: {'voc_name': 'Iso-propylbenzene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7406},
        43256: {'voc_name': 'Alpha-pinene', 'molecular_formula': 'C10H16', 'molecular_weight': 136.23, 'carbon_content': 10, 'hydrogen_content': 16, 'cid': 6654},
        45209: {'voc_name': 'N-propylbenzene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7668},
        45212: {'voc_name': 'M-ethyltoluene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 12100},
        45213: {'voc_name': 'P-ethyltoluene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 12160},
        45207: {'voc_name': '1,3,5-tri-m-benzene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7947},
        45211: {'voc_name': 'O-ethyltoluene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 11903},
        43257: {'voc_name': 'Beta-pinene', 'molecular_formula': 'C10H16', 'molecular_weight': 136.23, 'carbon_content': 10, 'hydrogen_content': 16, 'cid': 14896},
        45208: {'voc_name': '1,2,4-tri-m-benzene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 7247},
        43238: {'voc_name': 'N-decane', 'molecular_formula': 'C10H22', 'molecular_weight': 142.28, 'carbon_content': 10, 'hydrogen_content': 22, 'cid': 15600},
        45225: {'voc_name': '1,2,3-tri-m-benzene', 'molecular_formula': 'C9H12', 'molecular_weight': 120.19, 'carbon_content': 9, 'hydrogen_content': 12, 'cid': 10686},
        45218: {'voc_name': 'M-diethylbenzene', 'molecular_formula': 'C10H14', 'molecular_weight': 134.22, 'carbon_content': 10, 'hydrogen_content': 14, 'cid': 8864},
        45219: {'voc_name': 'P-diethylbenzene', 'molecular_formula': 'C10H14', 'molecular_weight': 134.22, 'carbon_content': 10, 'hydrogen_content': 14, 'cid': 7734},
        43954: {'voc_name': 'N-undecane', 'molecular_formula': 'C11H24', 'molecular_weight': 156.31, 'carbon_content': 11, 'hydrogen_content': 24, 'cid': 14257},
        43141: {'voc_name': 'N-dodecane', 'molecular_formula': 'C12H26', 'molecular_weight': 170.33, 'carbon_content': 12, 'hydrogen_content': 26, 'cid': 8182}
    }
    pass


class BaseSample:
    def __init__(self, front_path, back_path):
        self.front = FrontDetectorChromatogram(front_path)
        self.back = BackDetectorChromatogram(back_path)
# Specific sample types inherit from BaseSample
class Ambient(BaseSample):
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
        
