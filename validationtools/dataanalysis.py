# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 14:13:47 2025

@author: aengstrom
"""
#%% imports
import numpy as np
import pandas as pd
import pubchempy as pcp

#%% Reference Tables
#%%% Compound chemical categories
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
def sort_by_type(compound_categories):
    compound_df = pd.DataFrame(list(compound_categories.items()), columns=['compound', 'category'])
    group_by_category = compound_df.groupby(['category'])['compound'].apply(lambda x: [aqs_compound_to_code[c] for c in x])
    return group_by_category
#%%% MDL
rb_mdls = {
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
bv_mdls = {
    "Ethane": 0.0711,
    "Ethylene": 0.0642,
    "Propane": 0.0561,
    "Propylene": 0.0668,
    "Iso-butane": 0.0793,
    "N-butane": 0.0546,
    "Acetylene": 0.0512,
    "Trans-2-butene": 0.0409,
    "1-butene": 0.0604,
    "Cis-2-butene": 0.0543,
    "Cyclopentane": 0.0526,
    "Iso-pentane": 0.0583,
    "N-pentane": 0.0324,
    "1,3-butadiene": 0.0366,
    "Trans-2-pentene": 0.0484,
    "1-pentene": 0.032,
    "Cis-2-pentene": 0.0506,
    "2,2-dimethylbutane": 0.0529,
    "2,3-dimethylbutane": 0.0448,
    "2-methylpentane": 0.042,
    "3-methylpentane": 0.0486,
    "Isoprene": 0.0293,
    "2-methyl-1-pentene": 0.0512,
    "1-hexene": 0.0493,
    "N-hexane": 0.0593,
    "Methylcyclopentane": 0.0672,
    "2,4-dimethylpentane": 0.0848,
    "Benzene": 0.0525,
    "Cyclohexane": 0.0747,
    "2-methylhexane": 0.0896,
    "2,3-dimethylpentane": 0.0609,
    "3-methylhexane": 0.0586,
    "2,2,4-trimethylpentane": 0.0535,
    "N-heptane": 0.0615,
    "Methylcyclohexane": 0.0555,
    "2,3,4-trimethylpentane": 0.0478,
    "Toluene": 0.0674,
    "2-methylheptane": 0.0495,
    "3-methylheptane": 0.0778,
    "N-octane": 0.1171,
    "Ethylbenzene": 0.0619,
    "M&p-xylene": 0.2749,
    "Styrene": 0.0546,
    "O-xylene": 0.0586,
    "N-nonane": 0.0393,
    "Iso-propylbenzene": 0.0533,
    "Alpha-pinene": 0.2124,
    "N-propylbenzene": 0.0652,
    "M-ethyltoluene": 0.0953,
    "P-ethyltoluene": 0.0285,
    "1,3,5-tri-m-benzene": 0.0636,
    "O-ethyltoluene": 0.2738,
    "Beta-pinene": 0.1743,
    "1,2,4-tri-m-benzene": 0.2596,
    "N-decane": 0.0656,
    "1,2,3-tri-m-benzene": 0.067,
    "M-diethylbenzene": 0.077,
    "P-diethylbenzene": 0.1857,
    "N-undecane": 0.1558,
    "N-dodecane": 0.329
}
#%%% AQS upload file format
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
#%%% AQS compound codes
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
    "TNMHC": 43000
}

aqs_code_to_compound = {43202: 'Ethane', 
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
                               43000: 'TNMHC'}
def reverse_aqs_compound_to_code(aqs_code_dict = {
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
#%%% Compounds by column
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


#%% AQS txt file parsers
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
    aqs_code_to_voc_name = reverse_aqs_compound_to_code()
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
                aqs_code_to_voc_name = reverse_aqs_compound_to_code()
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

#%% Screening checks, based on reccomended screening checks in PAMS TAD Table 10-1, page 186
def check_totals(voc_df):
    required_cols = ['DATE/TIME', '43102', '43000']
    for col in required_cols:
        if col not in voc_df.columns:
            raise ValueError(f"Missing required column: {col}")

    totals_df = voc_df[required_cols]

    # Avoid division by zero
    nonzero = totals_df['43000'] != 0

    conditions = [
        (totals_df['43000'] <= 5, 'tnmhc_lt_5'),
        (totals_df['43000'] >= 500, 'tnmhc_gt_500'),
        (totals_df['43000'] >= totals_df['43102'], 'tnmhc_gt_tnmtc'),
        (((totals_df['43000'] - totals_df['43102']) / totals_df['43000'] >= 0.15) & nonzero, 'tnmtc_tnmhc_diff_gt_15'),
    ]

    flagged = []
    for cond, label in conditions:
        subset = totals_df.loc[cond].copy()
        subset['screen_reason'] = label
        subset['compounds'] = 'totals'
        flagged.append(subset)

    totals_check = pd.concat(flagged, ignore_index=True)
    return totals_check[['DATE/TIME', 'screen_reason','compounds']].copy()
def check_outliers(voc_df, summary_stats_df):
    # convert compound ID number column names to compound names
    voc_df_names = voc_df.copy()
    voc_df_names.columns = [voc_df_names.columns[0]] + [aqs_code_to_compound[int(c)].upper() for c in voc_df_names.columns[1:]]
    voc_df_names.iloc[:,1:] = voc_df_names.iloc[:,1:].apply(pd.to_numeric)
    summary_stats_df.rename(columns = {'Date':'statistic'},inplace=True)
    stats_clean = summary_stats_df.loc[summary_stats_df['statistic'].isin(["Mean", "Deviation"])].copy()
    stats_clean.iloc[:,1:] = stats_clean.iloc[:,1:].apply(pd.to_numeric)
    upper = stats_clean.loc[stats_clean["statistic"] == "Mean"].iloc[0, 1:] + 3 * stats_clean.loc[stats_clean["statistic"] == "Deviation"].iloc[0, 1:]
    lower = stats_clean.loc[stats_clean["statistic"] == "Mean"].iloc[0, 1:] - 3 * stats_clean.loc[stats_clean["statistic"] == "Deviation"].iloc[0, 1:]

    # Create rows with labels
    upper_row = pd.DataFrame([["upper", *upper.values]], columns=stats_clean.columns)
    lower_row = pd.DataFrame([["lower", *lower.values]], columns=stats_clean.columns)

    # Append them
    stats_clean = pd.concat([stats_clean, upper_row, lower_row], ignore_index=True)
    upper = stats_clean.set_index("statistic").loc["upper"]
    lower = stats_clean.set_index("statistic").loc["lower"]
    

    mask = voc_df_names.drop(columns=['DATE/TIME', 'TNMHC', 'TNMTC']).gt(upper) | voc_df_names.drop(columns=['DATE/TIME', 'TNMHC', 'TNMTC']).lt(lower)
    flags = pd.concat([voc_df_names['DATE/TIME'], mask], axis=1)
    mask_any = flags.drop(columns="DATE/TIME").any(axis=1)
    flags_filtered = flags.loc[mask_any].copy()
    flags_filtered["compounds"] = flags_filtered.drop(columns="DATE/TIME").apply(
    lambda row: [compound.capitalize() for compound in list(row.index[row])], axis=1
    )
    summary = flags_filtered[["DATE/TIME", "compounds"]].copy()
    summary.insert(loc = 1, column = 'screen_reason', value = 'outliers')
    return summary
def check_ratios(voc_df, mdls):
    mdl_number = {aqs_compound_to_code[name]: mdl for name, mdl in mdls.items()}
    voc_df.iloc[:,1:] = voc_df.iloc[:,1:].apply(pd.to_numeric)
    compound_by_type = sort_by_type(compound_categories = compound_categories)
    alkanes = [str(compound) for compound in compound_by_type['ALKANE']]
    alkenes = [str(compound) for compound in compound_by_type['ALKENE']]
    conditions = [
        ((voc_df['45201'] > voc_df['45202']) & (voc_df['45201'] > 3*mdl_number[45201]), 'benzene_gt_toluene',[45201,45202]),
        ((voc_df['45201'] > voc_df['43202']) & (voc_df['45201'] > 3*mdl_number[45201]), 'benzene_gt_ethane', [45201,43202]),
        ((voc_df['43203'] > voc_df['43202']) & (voc_df['43203'] > 3*mdl_number[43203]), 'ethylene_gt_ethane', [43203,43202]),
        ((voc_df['43205'] > voc_df['43204']) & (voc_df['43205'] > 3*mdl_number[43205]), 'propylene_gt_propane', [43205,43204]),
        ((voc_df['45204'] > voc_df['45109']) & (voc_df['45204'] > 3*mdl_number[45204]), 'oxylene_gt_mpxylene', [45204,45109]),
        ((voc_df['43263'] < voc_df['43291']) & (voc_df['43291'] > 3*mdl_number[43291]), '23dimethylpentane_gt_2methylhexane', [43263, 43291]),
        ((voc_df['43262'] < voc_df['43247']) & (voc_df['43247'] > 3*mdl_number[43247]), '24dimethylpentane_gt_methylcyclopentane', [43262, 43247]),
        ((voc_df['43214'] > voc_df['43212']) & (voc_df['43214'] > 3*mdl_number[43214]), 'isobutane_gt_nbutane', [43214, 43212]),
        ((voc_df['43230'] > .6*voc_df['43285']) & (voc_df['43230'] > 3*mdl_number[43230]), '3methylpentane_gt_2methylpentane', [43230,43285]),
        ((voc_df['43954'] > voc_df['43141']) & (voc_df['43954'] > 3*mdl_number[43954]), 'nundecane_gt_ndecane', [43954, 43141]),
        (~((voc_df['43221'] > voc_df['43220']) & (voc_df['43220'] > voc_df['43242'])) & (voc_df['43221'] > 3*mdl_number[43221]) & (voc_df['43220'] > 3*mdl_number[43220]) & (voc_df['43242'] > 3*mdl_number[43242]), 
         'not_isopentane_gt_npentane_gt_cyclopentane', [43221,43220,43242]),
        (voc_df[alkenes].sum(axis = 1) > voc_df[alkanes].sum(axis = 1), 'alkenes_gt_alkanes', ['alkanes', 'alkenes']) 
         ]
    flagged = []
    for cond, label, compounds in conditions:
        subset = voc_df.loc[cond].copy()
        subset['screen_reason'] = label
        subset['compounds'] = [[aqs_code_to_compound[compound] for compound in compounds if compound in aqs_code_to_compound.keys()]]*len(subset)
        flagged.append(subset)

    ratios_check = pd.concat(flagged, ignore_index=True)
    ratios = ratios_check[['DATE/TIME', 'screen_reason','compounds']].copy()
    return ratios

def screening_check(ambient_csv, summary_stats_cvs, mdls):
    summary_stats_df = pd.read_csv(summary_stats_cvs, header = 0, nrows = 5)
    voc_df = pd.read_csv(ambient_csv, header =2, parse_dates = [0])
    totals = check_totals(voc_df)
    outliers = check_outliers(voc_df, summary_stats_df)
    ratios = check_ratios(voc_df, mdls)
    all_checks = [totals, outliers, ratios]
    all_checks_df = pd.concat(all_checks, ignore_index=True)
    checks_by_date = all_checks_df.groupby('DATE/TIME').agg(list)
    return ratios
    
#%% Main
if __name__ == "__main__":
    compound_by_type = sort_by_type(compound_categories = compound_categories)
    test = screening_check(ambient_csv = r"D:\BV\working\validation\preprocessed_data\amount_crosstab_run_[S] - Copy.csv",
                                    summary_stats_cvs = r"D:\BV\working\validation\bv_summary_stats_2024.csv",
                                    mdls = bv_mdls)


    

        

    
    
    





