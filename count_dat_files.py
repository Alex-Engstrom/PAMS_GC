# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 10:15:23 2025

@author: aengstrom
"""
import os
from pathlib import Path
import re
from datetime import datetime, timedelta 
import shutil
def hours_since_year_start(target_date=None):
    """
    Calculate hours from start of year to target_date
    If target_date is None, uses current datetime
    """
    if target_date is None:
        target_date = datetime.now()
    
    year_start = datetime(target_date.year, 1, 1)
    time_difference = target_date - year_start
    total_hours = time_difference.total_seconds() / 3600
    
    return total_hours
def filename_info_extractor(filename):
    filename_pattern = r"(?P<site>[A-Z]{2})(?P<sample_type>[A-Z])(?P<month>[A-Z])(?P<day>\d{2})(?P<hour>[A-Z])"
    filename_match = re.match(filename_pattern, filename, re.IGNORECASE)
    if filename_match:
        return filename_match.groupdict()
    else:
        raise ValueError("Pattern not found in the filename.")
def verify_all_dat_files(input_folder):
    ''' function that searches through an unzipped gc daily data folder, and checks that the full day of .dat files are present
    inputs
    input_folder: path to daily gc data folder'''
    return_dict = {'unmatched_files' : [],
                   'missing_hours' : [],
                   'duplicate_hours' : [],
                   'incorrect_hours' : []}
    folder_path = Path(input_folder)
    #create list of files ending with .dat
    files = [file_path for file_path in folder_path.glob('*.dat')]
    #sorts files by last letetr of the filename stem (no .dat), should sort alphabetically from a to x for full day of data
    sorted_files = sorted(files, key = lambda s: s.stem[-1])
    #determines the year, month, day, and site from the foldername
    folder_pattern = r"([A-Z])(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})(?P<site>[A-Z]{2})"
    folder_match = re.match(folder_pattern, folder_path.stem, re.IGNORECASE)
    if folder_match:
        folder_dict = folder_match.groupdict()
    else:
        raise ValueError("Pattern not found in the folder name.")
    #dictionary to convert Orsat naming convention to 2 digit monthly naming conventions
    letter_to_number = {'a': '01', 'b': '02', 'c': '03', 'd': '04', 'e': '05', 
 'f': '06', 'g': '07', 'h': '08', 'i': '09', 'j': '10',
 'k': '11', 'l': '12', 'm': '13', 'n': '14', 'o': '15',
 'p': '16', 'q': '17', 'r': '18', 's': '19', 't': '20',
 'u': '21', 'v': '22', 'w': '23', 'x': '24', 'y': '25', 
 'z': '26'}
    hourly_files_check = {'a': [], 'b': [], 'c': [], 'd': [], 'e': [], 'f': [], 
 'g': [], 'h': [], 'i': [], 'j': [], 'k': [], 'l': [], 
 'm': [], 'n': [], 'o': [], 'p': [], 'q': [], 'r': [], 
 's': [], 't': [], 'u': [], 'v': [], 'w': [], 'x': []}
    for filename in sorted_files:
        if filename_info_extractor(filename.stem):
            filename_dict = filename_info_extractor(filename.stem)
            filename_month_digit = letter_to_number[filename_dict['month'].lower()]
            if filename_month_digit != folder_dict['month']:
                if filename not in return_dict['incorrect_hours']:
                    return_dict['incorrect_hours'].append(filename)
                
            #check that day on file matches day on folder
            if filename_dict['day'] != folder_dict['day']:
                if filename not in return_dict['incorrect_hours']:
                    return_dict['incorrect_hours'].append(filename)
            #get the hour value from the filename and make lowercase to match the hourly_files_check format
            hour = filename_dict['hour'].lower()
            if hour in hourly_files_check.keys():
                if filename_dict['day'] == folder_dict['day']:
                    hourly_files_check[hour].append(filename)               
            else:
                if filename not in return_dict['incorrect_hours']:
                    return_dict['incorrect_hours'].append(filename)
        else:
            return_dict['unmatched_files'].append(filename)
                
        #Iterates through hourly_files_check to check for missing or duplicate hours
        # Iterates through hourly_files_check to check for missing or duplicate hours
    for key, item in hourly_files_check.items():
        if not item:  # Check if the list is empty
            return_dict['missing_hours'].append(key)
        elif len(item) > 1:  # Check for duplicate hours
            return_dict['duplicate_hours'].append((key, item))
    return return_dict, hourly_files_check

  
def filter_and_sort_dat_files(input_folder, output_folder):
    letter_to_number = {'a': '00', 'b': '01', 'c': '02', 'd': '03', 'e': '04', 
 'f': '05', 'g': '06', 'h': '07', 'i': '08', 'j': '09',
 'k': '10', 'l': '11', 'm': '12', 'n': '13', 'o': '14',
 'p': '15', 'q': '16', 'r': '17', 's': '18', 't': '19',
 'u': '20', 'v': '21', 'w': '22', 'x': '23', 'y': '24', 
 'z': '25'}
    
    # Gets dictionary of files that dont belong, and dictionary with each hour of data
    bad_files, good_files = verify_all_dat_files(input_folder)
    
    for key, value in good_files.items():
        # Check if there are any files for this hour
        if value:  # This checks if the list is not empty
            fullpath = value[0]
            filename = fullpath.stem
            filename_info_dict = filename_info_extractor(filename)
            month_letter = filename_info_dict['month'].lower()
            month = letter_to_number[month_letter]
            day = filename_info_dict['day']
            hour_letter = filename_info_dict['hour'].lower()
            hour = letter_to_number[hour_letter]
            date = datetime(year=2025, month=int(month), day=int(day), hour=int(hour))
            tag = hours_since_year_start(target_date=date)
            print(tag)
            set_file_tags_using_timestamps(value[0], tag)
            shutil.copy2(fullpath, output_folder)
            
        else:
            print(f"No files found for hour {key}")
        
        
    
    


def set_file_tags_using_timestamps(file_path, tag):
    """Use file timestamps to encode sorting information"""
    
    tag_int = int(tag) if isinstance(tag, float) else tag
    
    try:
        # Use the tag value to set a specific modification time
        # This creates a predictable sorting order based on timestamps
        base_time = datetime(2025, 1, 1)  # Fixed base date
        target_time = base_time + timedelta(hours=tag_int)
        
        # Convert to timestamp
        timestamp = target_time.timestamp()
        
        # Set file modification time
        os.utime(file_path, (timestamp, timestamp))
        
        print(f"Set timestamp for {file_path.stem}: hour_{tag_int}")
        return True
        
    except Exception as e:
        print(f"Timestamp method failed for {file_path.stem}: {e}")
        return False

        
        


        
    

if __name__ == "__main__":
    main_folder_path = r"C:\Users\aengstrom\AutoGC\RB\validation\08\raw\unzipped"
    subfolders = os.listdir(main_folder_path)
    print(subfolders)
    month_data_dict = {}
    for subfolder in subfolders:
        # print(subfolder)
        # bad, good = verify_all_dat_files(os.path.join(main_folder_path, subfolder))
        # print(bad)
        filter_and_sort_dat_files(os.path.join(main_folder_path, subfolder), r"C:\Users\aengstrom\AutoGC\RB\validation\08\processed\dat_files")

    
