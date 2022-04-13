#!/usr/bin/env python3
"""
Count barcodes from a TRE-MPRA fastq.gz file.

(This is essentially a less robust version of starcode.)
"""
from typing import Dict, List
import io
import sys
import itertools
from collections import Counter
import gzip
import logging
import argparse

logging.getLogger().setLevel(logging.DEBUG)

# TODO: make this a path
def open_fastq(fname: str) -> io.TextIOWrapper:
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
        

def count_barcodes(fastqs: List[io.TextIOWrapper], barcode_length: int) -> Dict[str, int]:
    barcodes = []
    for fastq in fastqs:
        sequences = itertools.islice(fastq, 1, None, 4)
        barcodes.append(map(lambda x: x[:barcode_length], sequences))
    barcodes = itertools.chain.from_iterable(barcodes)
    return Counter(barcodes)


def main(seq_files, output, **kwargs):
    fastqs = [open_fastq(fname) for fname in seq_files]
    barcode_counts = count_barcodes(fastqs, 24)
    logging.info(f"Writing output to file '{output}'...")
    with gzip.open(output, 'wt') as output:
        for bc, n in barcode_counts.items():
           output.write(f"{bc}\t{n}\n") 
    logging.info(f"...done.")
    
parser = argparse.ArgumentParser(description="Count barcodes.")
parser.add_argument("seq_files", metavar="seq_file", nargs="+",
                    help="paths to demultiplexed FASTQ files to count from")
parser.add_argument("-o", "--output", metavar="FILE", required=True,
                    help="Provide a filename for output files.")

if __name__ == "__main__":
    args = vars(parser.parse_args())
    from pprint import pprint
    pprint(args)
    x = main(**args)

