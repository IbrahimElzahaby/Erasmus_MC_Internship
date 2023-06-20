#!/usr/bin/env python3
"""
Script other: Ibrahim ElZahaby
Script function: merge tsv files (plasmidfinder output) in one single file
Usage: python3 cat_tsv.py
"""

# import statement
import os

# declare output file variable
output_file = "all_TSVs.tsv"

# allocate the files path
directory = "/mnt/DATAPOOL/mmibstudentnew/DATA/plasmidfinder_output/tsv_files"
header_written = False

# open tsv files
with open(output_file, "w") as outfile:
    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):
            file_path = os.path.join(directory, filename)

            with open(file_path, "r") as infile:
                lines = infile.readlines()
                # Check if the file has more than one line (excluding the header)
                if len(lines) > 1:
                    if not header_written:
                        # Write the header line only once
                        outfile.write(lines[0])
                        header_written = True
                    # Write the remaining lines excluding the header
                    outfile.writelines(lines[1:])

