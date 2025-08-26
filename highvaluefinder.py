# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 13:34:00 2025

@author: aengstrom
"""

import pandas as pd
def open_and_sort_sample_data(data_path, upper_calibration_point):
    """
    Find all instances where values exceed the threshold.
    
    Parameters:
    data_path: CSV from Orsat Crosstab CSV function
    upper_calbration_point: upper point in calibration curve (ppbC), searches for values above this threshhold
    
    Returns:
    list of tuples: (column_name, datetime, value)
    """
    # Read the CSV without headers first
    df = pd.read_csv(data_path, header=None)
    
    # Extract the specific headers
    first_col_header = df.iloc[2, 0]  # Row 3, Column 1 (0-indexed: row 2, column 0)
    other_col_headers = df.iloc[1, 1:].tolist()  # Row 2, Columns 2+ (0-indexed: row 1, columns 1+)
    
    # Combine the headers
    column_names = [first_col_header] + other_col_headers
    
    # Re-read the CSV skipping the first 3 rows and use custom headers
    df = pd.read_csv(data_path, header=None, skiprows=3)
    df.columns = column_names

    exceedances = []  # Initialize empty list to store results
    
    # Get the datetime column (assuming it's the first column)
    datetime_col = df.columns[0]  # Access the first column name
    
    # Iterate through all columns except the first one
    for column in df.columns[1:-2]:  # Loop starts from second column
        # Create boolean mask for values exceeding threshold
        mask = df[column] > upper_calibration_point  # Returns True/False Series
        
        # Get the rows where condition is True
        exceeding_rows = df.loc[mask]  # Filter DataFrame using boolean mask
        
        # Append results to the list
        for index, row in exceeding_rows.iterrows():
            exceedances.append((column, row[datetime_col], row[column]))
    
    return exceedances
if __name__ == "__main__":
    exceedance= open_and_sort_sample_data(r"C:\Users\aengstrom\Downloads\amount_crosstab_run_[S].csv", upper_calibration_point=45)
    for j in exceedance:
        print(j)
    