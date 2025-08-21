# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 12:02:26 2025

@author: aengstrom
"""

import os
import calendar

def create_date_folders():
    # Ask the user for the root directory
    root_dir = input("Enter the root directory where you want the folders created: ").strip('"')

    # Ask the user for the year
    year = int(input("Enter the year to create the folder structure for: "))

    # Create the root year folder
    year_path = os.path.join(root_dir, str(year))
    os.makedirs(year_path, exist_ok=True)
    print(f"Created folder: {year_path}")

    # Loop through all 12 months
    for month in range(1, 13):
        month_path = os.path.join(year_path, f"{month:02d}")
        os.makedirs(month_path, exist_ok=True)
        print(f"  Created month folder: {month_path}")

        # Get correct number of days in this month
        num_days = calendar.monthrange(year, month)[1]

        # Create day folders
        for day in range(1, num_days + 1):
            day_path = os.path.join(month_path, f"{day:02d}")
            os.makedirs(day_path, exist_ok=True)

        print(f"    Created {num_days} day folders for {month:02d}/{year}")

    print("\n Folder structure created successfully!")

# Run the function
if __name__ == "__main__":
    create_date_folders()