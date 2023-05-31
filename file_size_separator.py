#!/usr/bin/env python3
"""
Script other: Ibrahim ElZahaby
Student number: 1069624
Script function: separate 239 gff files above 1 kb and below 1 kb into two folders
Usage: python3 file_size_separator.py
"""
import os
import shutil

# Define the paths to the input and output directories
gffs_dir = '/mnt/DATAPOOL/mmibstudentnew/dummy/output_files/panaroo_data/'
large_gffs = '/mnt/DATAPOOL/mmibstudentnew/dummy/output_files/panaroo_data/large_files/'
small_gffs = '/mnt/DATAPOOL/mmibstudentnew/dummy/output_files/panaroo_data/small_files/'

# Loop over the files in the input directory
for filename in os.listdir(gffs_dir):
    # Construct the full path to the file
    filepath = os.path.join(gffs_dir, filename)
    # Check if the file size is greater than 5 KB
    if os.path.getsize(filepath) > 1 * 1024:
        # Move the file to the large files directory
        shutil.move(filepath, large_gffs)
    else:
        # Move the file to the small files directory
        shutil.move(filepath, small_gffs)
