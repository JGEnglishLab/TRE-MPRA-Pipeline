#!/usr/bin/env python3
"""
Map barcodes from the TRE-MPRA library architecture:

    R2 primer
    ________\\
        MluI     oligo    XbaI   barcode
    ---ACGCGT|---------/-|TCTAGA|BDHVBDHVBDHVHVDBHVDBHVDB|
    ---TGCGCA|---------/-|AGATCT|BDHVBDHVBDHVHVDBHVDBHVDB|
                                                          ___________
                                                          \ R1 primer

Warning: the Levenshtein distance search is really inefficient here.
"""
import sys
import re
import argparse
import logging
import itertools
from collections import Counter, namedtuple, defaultdict
import gzip
from Bio import SeqIO
import numpy as np
import Levenshtein

LibMember = namedtuple("LibMember", ['id', 'distance', 'bc_counts'])

logging.getLogger().setLevel(logging.DEBUG)

# TODO: make this a path
def open_fastq(fname: str): 
    try:
        fastq = gzip.open(fname, "rt")
        # try to read one byte
        fastq.read(1) 
        logging.debug(f"{fname} is gzip file")
        fastq.seek(0)
    except gzip.BadGzipFile:
        logging.debug(f"{fname} is not gzip file")
        fastq = open(fname)
    return fastq

# core oligo sequence flanked by MluI
OLIGO_RE = re.compile(r"(.*)ACGCGT", flags=re.IGNORECASE)

def extract_oligo(rec):
    """select the relevant sequence between the MluI and XbaI sites"""
    # revcomp to put seq in orientation of R1 (see module docstring)
    seq = str(rec.seq.reverse_complement())
    match = OLIGO_RE.search(seq)
    if not(match):
        raise ValueError(f"Could not find oligo in {rec.id}")
    return match.group(1)

def make_reference(ref_name):
    library = { 
        extract_oligo(rec) : LibMember(rec.id, 0, Counter())
        for rec in SeqIO.parse(ref_name, "fasta")
    }
    return library

def find_oligo_barcode(seq):

    xba_i = seq.find("TCTAGA")
    return(seq[:xba_i], seq[xba_i+6:])
    if xba_i < 22 or xba_i > 24:
        return(None)
    else:
        pass

def map_barcodes(library, seq_files):
    orphans = defaultdict(Counter)
    for fname in seq_files:
        fastq = open_fastq(fname)
        # for rec in SeqIO.parse(fname, format="fastq"):
        for seq in itertools.islice(fastq, 1, None, 4):
            seq = seq.strip()
            # extract the oligo and barcode from the seq 
            extracted = find_oligo_barcode(seq)
            if extracted is None:
                continue
            bc, oligo = extracted
            
            match = library.get(oligo)
            if match:
                match.bc_counts[bc] += 1
            else:
                orphans[oligo][bc] += 1

    return(library, orphans)

def rescue_orphans(library, orphans, max_dist):
    mutants = {}
    unrescued = {}

    for orphan_oligo, barcodes in orphans.items():
        # a few heuristics
        if (orphan_oligo.count("N") + orphan_oligo.count("n") > max_dist
            or len(orphan_oligo) < 134 - max_dist):
            unrescued[orphan_oligo] = barcodes
            continue 
        
        min_dist = max_dist + 1
        nearest_oligo = None
        clash = False
        for lib_oligo in library:
            if abs(len(lib_oligo) - len(orphan_oligo)) > max_dist:
                continue
            dist = Levenshtein.distance(orphan_oligo, lib_oligo)
            if dist < min_dist:
                min_dist = dist
                nearest_oligo = lib_oligo
            elif dist == min_dist and nearest_oligo is not None:
                # if two members are of the same Levenshtein distance away, 
                # we treat this as a clash and remove the barcode
                clash = True
                
        if not clash and nearest_oligo is not None:
            mutants[orphan_oligo] = LibMember(library[nearest_oligo].id, min_dist, barcodes)
        else:
            unrescued[orphan_oligo] = barcodes

    return(mutants, unrescued)

def export_tsvs(library, orphans, prefix):
    with open(f"{prefix}_orphans.tsv", 'w') as outf:
        outf.write("sequence\tbarcode\toccurrence\n")
        for sequence, barcodes in orphans.items():
            for barcode, occurrence in barcodes.items():
                outf.write(f"{sequence}\t{barcode}\t{occurrence}\n")
    with open(f"{prefix}_mapped.tsv", 'w') as outf:
        outf.write("id\tdistance\tbarcode\toccurrence\n")
        for (_, (name, distance, barcodes)) in library.items():
            # TODO: add class and split fields
            for barcode, occurrence in barcodes.items():
                outf.write(f"{name}\t{distance}\t{barcode}\t{occurrence}\n")


def main(reference, seq_files, distance=0, output=None): 
    import os
    for fname in seq_files:
        assert os.path.exists(fname) 
    if output is None:
        output = "barcode" 

    logging.info(f"Reading reference file '{reference}'...")
    library = make_reference(reference)
    logging.info(f"Finding exact oligo matches in files {seq_files!r}...")
    library, orphans = map_barcodes(library, seq_files)
    
    if distance > 0:
        logging.info(f"Recovering orphan barcodes within a distance of {distance}...")
        rescued, orphans = rescue_orphans(library, orphans, max_dist=distance)
        library.update(rescued)
     
    export_tsvs(library, orphans, output)  
    return library, orphans
     

parser = argparse.ArgumentParser(description="Map barcodes.")
parser.add_argument("reference", 
                    help="Path to FASTA-formatted library reference")
parser.add_argument("seq_files", metavar="seq_file", nargs="+",
                    help="paths to demultiplexed FASTQ files to map from")
parser.add_argument("-d", "--distance",  metavar="N", type=int, default=0,
                    help="Retain barcodes whose oligo sequence are an edit distance <= N away from reference")
parser.add_argument("-o", "--output", metavar="PREFIX",
                    help="Provide a prefix for output files. {PREFIX}_orphans.tsv, {PREFIX}_mapped.tsv, etc.")


if __name__ == "__main__":
    args = vars(parser.parse_args())
    from pprint import pprint
    pprint(args)
    x = main(**args)

