#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-13 17:11:58
@LastEditTime : 2020-05-20 15:09:21
@Description  : adjust read strandedness
                only output reads contain 3'adapter
'''


import pandas as pd
import pysam
import click
from Bio import SeqIO
import gzip


@click.command()
@click.option('-i', '--infile', required=True, help='Input fastq file.', type=click.Path(exists=True))
@click.option('-p', '--polya_file', required=True, help='Input polya_file file.', type=click.Path(exists=True))
@click.option('-o', '--outfile', required=True, help='Output adjust file.')
def main(infile, outfile, polya_file):
    # load tail_length and read_type data
    df = pd.read_csv(polya_file, sep='\t', index_col=0)
    polya_type = df['read_type'].map(lambda x: False if x=='polyT' else True).to_dict()
    
    with gzip.open(infile, 'rt') as f, gzip.open(outfile, 'wt') as o:
        for read in SeqIO.parse(f, 'fastq'):
            if read.id in polya_type:
                if not polya_type[read.id]:
                    read.seq = read.seq.reverse_complement()
                SeqIO.write(read, o, 'fastq')
                
                
if __name__ == "__main__":
    main()