#!/usr/bin/env python3
"""
Script other: Ibrahim ElZahaby
Script function: retrieve remaining fasta files for the whole plasmid (RefSeq) counts from NCBI nucleotide database
Usage: python3 progress_fasta.py
"""

# Import statement
from Bio import Entrez, SeqIO

def fetch_sequences(accession_list):
    """
    Fetches the nucleotide sequences for plasmids RefSeq accession list from the NCBI nucleotide database.
    Parameters:
        accession_list (str): The accession numbers of the sequences.
    Returns:
        records: The fasta records containing the fetched sequences.
    """
    # Set email address
    Entrez.email = "ibrahimgolden342@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=accession_list, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    return records

# Read accession numbers from the sequence.seq file
acc_names = "sequence.seq"

with open(acc_names, "r") as file:
    accession_list = [line.strip() for line in file]

# load progress files information
progress_file = "progress.txt"
processed_accessions = []
try:
    with open(progress_file, "r") as progress:
        processed_accessions = [line.strip() for line in progress]
except FileNotFoundError:
    pass

# Fetch sequences and save in separate fasta files
for accession in accession_list:
    if accession not in processed_accessions:
        try:
            sequence_record = fetch_sequences(accession)
            output_filename = f"{accession}.fasta"
            with open(output_filename, "w") as output_file:
                SeqIO.write(sequence_record, output_file, "fasta")
            processed_accessions.append(accession)
        except Exception as e:
            print(f"Error fetching sequence for accession {accession}: {e}")
    
# Save progress
with open(progress_file, "w") as progress:
    for accession in processed_accessions:
        progress.write(accession + "\n")
