#!/usr/bin/env python3
"""
Script other: Ibrahim ElZahaby
Script function: retrieve fasta files for the whole plasmid (RefSeq) counts from NCBI nucleotide database
Usage: python3 retrieve_fasta.py
"""

# Import statement
from Bio import Entrez, SeqIO

def fetch_sequences(accession_list):
    """
    Fetches the nucleotide sequences for a given accession number from the NCBI nucleotide database.
    Parameters:
        accession_list (str): The accession number of the sequences.
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

# Fetch sequences and save in separate fasta files
for accession in accession_list:
    fasta_records = fetch_sequences(accession)
    output_fasta = f"{accession}.fasta"
    with open(output_fasta, "w") as output_fasta:
        SeqIO.write(fasta_records, output_fasta, "fasta")


