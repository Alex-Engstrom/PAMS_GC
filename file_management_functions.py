# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 13:33:27 2025

@author: aengstrom
"""

import os
import calendar
import re
import glob

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

import os
import zipfile
import shutil
from pathlib import Path

def unzip_files(source_directory, destination_directory, delete_zip_after_extract=False, create_subfolders=True):
    """
    Unzips all zip files in a source directory to a destination directory.
    
    Parameters:
    source_directory (str): Path to the directory containing zip files
    destination_directory (str): Path where extracted files should be saved
    delete_zip_after_extract (bool): Whether to delete the zip file after extraction (default: False)
    create_subfolders (bool): Whether to create subfolders for each zip file (default: True)
    
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
                # Create subfolder with zip file name (without .zip extension)
                extract_folder_name = zip_file.stem  # removes .zip extension
                extract_path = dest_dir / extract_folder_name
            else:
                # Extract directly to destination directory
                extract_path = dest_dir
            
            # Create extraction directory if it doesn't exist
            extract_path.mkdir(parents=True, exist_ok=True)
            print(f"Extracting to: {extract_path}")
            
            # Extract the zip file
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                # Get list of files in zip
                file_list = zip_ref.namelist()
                print(f"Found {len(file_list)} file(s) in zip")
                
                # Extract all files
                zip_ref.extractall(extract_path)
                print(f"✓ Successfully extracted {len(file_list)} file(s)")
                
                # List extracted files
                for extracted_file in file_list:
                    print(f"  - {extracted_file}")
            
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


    
def move_dat_files(month_dir):
    file_dir = os.path.join(month_dir, 'unzipped')
    files_in_dir = os.listdir(file_dir)
    sorted_files = sorted(files_in_dir, key = lambda s: s[-4:-2])
    print(sorted_files)
    week_folder_names = ['']
    dat_file_paths = glob.glob(file_dir+r"/**/*.dat", recursive = True)
    print(dat_file_paths)
    output_folder_path = os.path.join(month_dir, 'dat_file_dump')
    os.makedirs(output_folder_path, mode = 0o755, exist_ok = True)
    for path in dat_file_paths:
        shutil.copy2(path,output_folder_path)
        
        
   
def transfer_raw_gc_data(source_path = None, output_dir = r"C:\AutoGCData", station_initials= None ,year = None, month = None):
    output_validation_dir = os.path.join(output_dir, station_initials, 'validation')
    if source_path == None:
        base_dir = r"\\UTAH\DEQ\DAQ\SHARED\PLAN\AMC\GC-Agilent+Markes\Data\ORIGINAL"
        sites= os.listdir(base_dir)
        for sitename in sites:
            if station_initials != None:
                if station_initials in sitename:
                    source_path = os.path.join(base_dir, sitename, year, f"{year}{month}")
            else:
                print("Please ensure station initials, year, and month are entered properly")
    all_entries = os.listdir(source_path)
    all_entries_dicts = []
    pattern = r"([A-Z])(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})(?P<site>[A-Z]{2}).(?P<filetype>zip)"
    for filename in all_entries:
        match = re.match(pattern, filename, re.IGNORECASE)
        if match:
            file_dict = match.groupdict()
            file_dict['filename'] = filename
            all_entries_dicts.append(file_dict)
    for file in all_entries_dicts:
        for folder in os.listdir(output_validation_dir):
            if re.fullmatch("\d{4}", folder):
                yr = folder
            for subfolder in os.listdir(output_validation_dir+f"\{yr}"):
                if re.fullmatch("\d{2}", subfolder):
                    mnth = subfolder
                    if '20'+file['year'] == yr and file['month'] == mnth:
                        try:
                            source_path_file = os.path.join(source_path, file['filename'])
                            if os.path.exists(source_path_file):
                                print('exists')
                            destination_dir = os.path.join(output_validation_dir, yr, mnth, 'zipped')
                            if os.path.exists(destination_dir) == False:
                                os.mkdir(destination_dir)
                            destination_path = os.path.join(destination_dir, file['filename'])
                            shutil.copyfile(source_path_file, destination_path)
                            print(f"File '{source_path_file}' copied to '{destination_path}' successfully.")
                        except FileNotFoundError:
                            print("Source file not found.")
                        except IsADirectoryError:
                            print("Destination is a directory, not a file.")
                        except PermissionError:
                            print("Permission denied to write to the destination.")
                        except shutil.SameFileError:
                            print("Source and destination are the same file.")
                        
                else:
                    continue
            else:
                continue
    source_directory = os.path.join(destination_dir)
    destination_directory = os.path.join(output_validation_dir, year, month, 'unzipped')
    unzip_files(source_directory, destination_directory, delete_zip_after_extract=False, create_subfolders=True)
            


if __name__ == "__main__":
    move_dat_files(month_dir=r"C:\AutoGCData\RB\validation\2025\08")


