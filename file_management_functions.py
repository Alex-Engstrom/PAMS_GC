# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 13:33:27 2025

@author: aengstrom
"""

import os
import calendar
import re
import pandas as pd
from pathlib import Path
#%% File structure generators
def create_date_folders(root_dir, site, subfolders,year, make_month = True, make_day = True):

    # Create the root year folder
    site_path = os.path.join(root_dir, site)
    os.makedirs(site_path, exist_ok=True)
    print(f"Created folder: {site_path}")
    for i, subfolder in enumerate(subfolders):
        subfolder_path = os.path.join(site_path, subfolder, year)
        os.makedirs(subfolder_path, mode = 0o755, exist_ok = True)
        if make_month[i] == True:
            # Loop through all 12 months
            for month in range(1, 13):
                month_path = os.path.join(subfolder_path, f"{month:02d}")
                os.makedirs(month_path, mode = 0o755, exist_ok=True)
                print(f"  Created month folder: {month_path}")
        
                # Get correct number of days in this month
                num_days = calendar.monthrange(int(year), month)[1]
                if make_day[i] == True:
                # Create day folders
                    for day in range(1, num_days + 1):
                        day_path = os.path.join(month_path, f"{day:02d}")
                        os.makedirs(day_path, mode = 0o774, exist_ok=True)
            
                    print(f"    Created {num_days} day folders for {month:02d}/{year}")

    print("\n Folder structure created successfully!")
def generate_monthly_final_folder_structure(root_dir, sitename, year, month):
    root = Path(root_dir)
    zipfile = root / f'{sitename}{str(year)}{month:02d}v1'
    zipfile.mkdir(exist_ok = True)
    folders = [zipfile / 'AQS',
               zipfile / 'FINAL',
               zipfile / 'FINAL' / 'dat_and_txt',
               zipfile / 'FINAL' / 'dat_and_txt' / 'Final',
               zipfile / 'FINAL' / 'dat_and_txt' / 'Original',
               zipfile / 'FINAL' / 'dat_and_txt' / 'RPO',
               zipfile / 'MDVR',
               zipfile / 'Original']
    for folder in folders:
        folder.mkdir(exist_ok = True)
import subprocess

        

    
        
    
#%% Unzipper Functions
import os
import zipfile
import shutil
from pathlib import Path

def unzip_files(source_directory, destination_directory, delete_zip_after_extract=False, create_subfolders=True):
    """
    Unzips all zip files in a source directory to a destination directory,
    preserving the original modification dates.
    
    Parameters:
    source_directory (str): Path to the directory containing zip files
    destination_directory (str): Path where extracted files should be saved
    delete_zip_after_extract (bool): Whether to delete the zip file after extraction
    create_subfolders (bool): Whether to create subfolders for each zip file
    
    Returns:
    list: List of successfully extracted zip files
    """
    
    # Convert to Path objects for easier handling
    source_dir = Path(source_directory)
    dest_dir = Path(destination_directory)
    
    # Check if source directory exists
    if not source_dir.exists():
        print(f"Error: Source directory does not exist: {source_directory}")
        return []
    
    if not source_dir.is_dir():
        print(f"Error: Source path is not a directory: {source_directory}")
        return []
    
    # Create destination directory if it doesn't exist
    dest_dir.mkdir(parents=True, exist_ok=True)
    print(f"Destination directory: {dest_dir}")
    
    # Find all zip files in source directory
    zip_files = list(source_dir.glob("*.zip"))
    print(f"Found {len(zip_files)} zip file(s) in source directory")
    
    successfully_extracted = []
    
    for zip_file in zip_files:
        try:
            print(f"\nProcessing: {zip_file.name}")
            
            # Determine extraction path
            if create_subfolders:
                extract_folder_name = zip_file.stem
                extract_path = dest_dir / extract_folder_name
            else:
                extract_path = dest_dir
            
            # Create extraction directory if it doesn't exist
            extract_path.mkdir(parents=True, exist_ok = True)
            print(f"Extracting to: {extract_path}")
            
            # Extract the zip file with date preservation
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                file_list = zip_ref.namelist()
                print(f"Found {len(file_list)} file(s) in zip")
                
                # Extract each file individually to preserve modification dates
                for file_info in zip_ref.infolist():
                    # Extract the file
                    zip_ref.extract(file_info, extract_path)
                    
                    # Get the full path to the extracted file
                    extracted_file_path = extract_path / file_info.filename
                    
                    # Preserve the original modification date
                    if not file_info.is_dir():
                        date_time = file_info.date_time
                        timestamp = datetime(*date_time).timestamp()
                        os.utime(extracted_file_path, (timestamp, timestamp))
                        
                    print(f"  - {file_info.filename} (mod date preserved)")
                
                print(f"✓ Successfully extracted {len(file_list)} file(s)")
            
            successfully_extracted.append(str(zip_file))
            
            # Delete zip file if requested
            if delete_zip_after_extract:
                zip_file.unlink()
                print(f"✓ Deleted original zip file: {zip_file.name}")
                
        except zipfile.BadZipFile:
            print(f"✗ Error: {zip_file.name} is not a valid zip file")
        except PermissionError:
            print(f"✗ Error: Permission denied when processing {zip_file.name}")
        except Exception as e:
            print(f"✗ Error processing {zip_file.name}: {e}")
    
    print(f"\n=== Summary ===")
    print(f"Total zip files processed: {len(zip_files)}")
    print(f"Successfully extracted: {len(successfully_extracted)}")
    print(f"Failed: {len(zip_files) - len(successfully_extracted)}")
    
    return successfully_extracted


   


#%% Move, extract, sort, and convert daily zip files

def transfer_raw_gc_data(source_dir, output_dir):
    '''
    Sorts through a folder of daily zip files, extracts and moves them to a new directory
    
    '''
    source_path = Path(source_dir)
    output_dir_path = Path(output_dir)
    all_entries = source_path.iterdir()
    all_entries_dicts = []
    pattern = r"([A-Z])(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})(?P<site>[A-Z]{2}).(?P<filetype>zip)"
    for filename in all_entries:
        match = re.match(pattern, filename.name, re.IGNORECASE)
        if match:
            file_dict = match.groupdict()
            file_dict['filename'] = filename
            all_entries_dicts.append(file_dict)
    unzipped_destination_dir = output_dir_path / 'unzipped'
    unzip_files(source_path, unzipped_destination_dir, delete_zip_after_extract=False, create_subfolders=True)
    return unzipped_destination_dir 

def move_dat_files(input_dir, output_dir):
    ''' 
    Sorts through folder of unzipped daily gc data folders, identifies the .dat files and moves them to a new folder.
    '''
    file_dir = Path(input_dir)
    files_in_dir = file_dir.iterdir()
    sorted_files = sorted(files_in_dir, key = lambda s: s.name[-4:-2])
    print(sorted_files)
    dat_file_paths = file_dir.rglob(r"*.dat")
    print(dat_file_paths)
    output_folder_path = Path(output_dir, 'dat_file_dump')
    output_folder_path.mkdir(mode = 0o755, exist_ok = True)
    for path in dat_file_paths:
        if path.stem not in [file.stem for file in list(output_folder_path.iterdir())]:
            shutil.copy2(path,output_folder_path)
        else:
            (output_folder_path / 'copies').mkdir(mode = 0o755, exist_ok = True)
            shutil.copy2(path,output_folder_path / 'copies')
            
    return output_folder_path
        

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
        filename = filename.stem
        if filename_info_extractor(filename):
            filename_dict = filename_info_extractor(filename)
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
  
def filter_and_sort_dat_files(input_folder, output_folder):
    letter_to_number = {'a': '00', 'b': '01', 'c': '02', 'd': '03', 'e': '04', 
 'f': '05', 'g': '06', 'h': '07', 'i': '08', 'j': '09',
 'k': '10', 'l': '11', 'm': '12', 'n': '13', 'o': '14',
 'p': '15', 'q': '16', 'r': '17', 's': '18', 't': '19',
 'u': '20', 'v': '21', 'w': '22', 'x': '23', 'y': '24', 
 'z': '25'}
    
    # Gets dictionary of files that dont belong, and dictionary with each hour of data
    input_folder_path = Path(input_folder)
    dat_files = input_folder_path.rglob('*.dat')
    
    for dat_file in dat_files:
        filename = dat_file.stem
        filename_info_dict = filename_info_extractor(filename)
        month_letter = filename_info_dict['month'].lower()
        month = letter_to_number[month_letter]
        day = filename_info_dict['day']
        hour_letter = filename_info_dict['hour'].lower()
        hour = letter_to_number[hour_letter]
        date = datetime(year=2025, month=int(month), day=int(day), hour=int(hour))
        tag = hours_since_year_start(target_date=date)
        print(tag)
        set_file_tags_using_timestamps(dat_file, tag)
        if input_folder != output_folder:
            shutil.copy2(dat_file, output_folder)
        else:
            continue
            

def convert_file_to_pdf(input_file_path, output_pdf_path):
    """
    Converts a DOCX file to PDF using LibreOffice
    """
    soffice_path = r"C:\Program Files\LibreOffice\program\soffice.exe"

    output_file_path = Path(output_pdf_path)
    try:
        temp_output_dir = output_pdf_path.parent
        subprocess.run(
            [
                soffice_path,
                "--headless",
                "--convert-to", "pdf",
                "--outdir", str(output_file_path.parent),
                str(input_file_path)
            ],
            check=True,
            capture_output=True,
            text=True
        )
        # Rename the converted file to the desired name
        temp_pdf_path = temp_output_dir / f"{input_file_path.stem}.pdf"
        
        if temp_pdf_path.exists():
            # Rename to the desired filename
            temp_pdf_path.rename(output_file_path)
            print(f"Successfully converted and renamed '{input_file_path}' to '{output_pdf_path}'.")
        else:
            print(f"Error: Converted file not found at expected location: {temp_pdf_path}")
        print(f"Successfully converted '{input_file_path}' to '{output_pdf_path}'.")
    except Exception as e:
        print(e)
        try: 
            temp_pdf_path.unlink()
        except OSError as e:
            print(f"Error deleting file '{temp_pdf_path}': {e}")
            



def convert_folder_contents_to_pdf(input_folder_dir, output_folder_dir):
    input_folder_path = Path(input_folder_dir)
    output_folder_path = Path(output_folder_dir)
    
    output_folder_path.mkdir(parents=True, exist_ok=True)
    extensions = ('*.docx', '*.xlsx', '*.xlsm')
    convertible_files = []
    for extension in extensions:
        convertible_files.extend(input_folder_path.rglob(extension))
    for file in convertible_files:
        file_stats = file.stat()
        modified_timestamp = file_stats.st_mtime
        modified_datetime = datetime.fromtimestamp(modified_timestamp)
        new_pdf_name = str(file.stem) + ' modified ' + modified_datetime.strftime("%Y-%m-%d-%H-%M-%S") + '.pdf'
        new_pdf_path = output_folder_path / new_pdf_name
        convert_file_to_pdf(file, output_pdf_path = new_pdf_path)
        
        
            
def set_up_validation_folder(source_dir, output_dir):
    unzipped_dir = transfer_raw_gc_data(source_dir = source_dir, output_dir = output_dir)
    
    print(str(unzipped_dir))
    dat_dump_dir = move_dat_files(input_dir = str(unzipped_dir), output_dir = output_dir)
    daily_summaries = []
    for daily_gc_folder in unzipped_dir.iterdir():
        daily_summaries.append(verify_all_dat_files(input_folder = daily_gc_folder))
    filter_and_sort_dat_files(input_folder = str(dat_dump_dir) , output_folder = str(dat_dump_dir))
    convert_folder_contents_to_pdf(input_folder_dir= unzipped_dir, output_folder_dir = output_dir)
    return daily_summaries
    
    
    
#%% File Naming Functions

def rename_dattxt_files_to_txt(input_directory, output_directory = None):
    # list the directory contents
    #directory = "U:\\PLAN\\AMC\\GC-Agilent+Markes\\Data\\MDVR\\2021\\202105 MDL\\"
    #directory = "U:\\PLAN\\Rkaullman\\GC Data Exp\\"
    #directory = "D:\\Testing\\"
    input_directory_path = Path(input_directory)
    if output_directory == None:
        output_directory_path = input_directory_path
    else:
        output_directory_path = Path(output_directory)
    if input_directory_path.exists():
        files = input_directory_path.iterdir()
    
        
        for filename in files:
            try:
                if filename.name.endswith(".dat.tx1"):
                    newname = filename.name.replace(".dat", "", 1)
                    print(filename.name + " => " + newname)
                    # Change hwqe20h.dattx1 to hwqe20h.tx1
                    filename.rename(output_directory_path / newname)
            except FileExistsError:
                print('file already exists')
    

#%% main
if __name__ == "__main__":
    aqs_df = generate_df_from_aqs_file(r"C:\Users\aengstrom\Desktop\RD_SITE-RB_RedButte_INSERT_08-01-2025_to_08-31-2025.txt")



