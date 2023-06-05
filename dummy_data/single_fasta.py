#!/usr/bin/env python3
"""
Script other: Ibrahim ElZahaby
Micro section number: 097220
Script function: python code generate 261 plasmid sequence records (output txt files) from fasta file
Usage: python3 single_fasta.py
"""
def split_fasta(fasta_file):
    """
    Splits a multi-sequence FASTA file into individual files, one for each sequence.
    The output files will be named with the accession number extracted from the sequence header.
    The accession number is assumed to be the first string after the third '|' character in the header.
    :param fasta_file: input FASTA file includes all genomic records
    """

# declare input fasta file variable
fasta_file = "ctx_seqs.fasta"
# open fasta file
with open(fasta_file, "r") as f:
    # split fasta file into individual sequences
    sequences = f.read().split(">")

    for seq in sequences[1:]:
        # split the sequence into lines
        lines = seq.strip().split("\n")
        # extract the header line
        header = lines[0]
        # join the remaining lines into a single sequence
        seq = "".join(lines[1:])
        # extract the accession number from the header
        accession = header.split("|")[3]

        # write the sequence to a text file with the accession number as the output file name
        with open(f"{accession}.txt", "w") as out_file:
            out_file.write(f">{header}\n{seq}")