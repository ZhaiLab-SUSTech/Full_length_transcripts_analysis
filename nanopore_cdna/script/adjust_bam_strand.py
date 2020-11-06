#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-07-13 10:09:35
@LastEditTime : 2020-07-13 10:50:38
@Description  : 
'''


import pysam
import click
import pandas as pd


@click.command()
@click.option('-i', '--infile', required=True, help='Input bam file.', type=click.Path(exists=True))
@click.option('-p', '--polya_file', required=True, help='Input polya_file file.', type=click.Path(exists=True))
@click.option('-o', '--outfile', required=True, help='Output adjust bam.')
def main(infile, outfile, polya_file):
    # load tail_length and read_type data
    df = pd.read_csv(polya_file, sep='\t', index_col=0)
    polya_type = df['read_type'].map(lambda x: False if x=='polyT' else True).to_dict()

    inbam = pysam.AlignmentFile(infile, 'rb')
    outbam = pysam.AlignmentFile(outfile, 'wb', template=inbam)
        
    for read in inbam:
        if read.query_name in polya_type:
            if not polya_type[read.query_name]:
                if read.is_reverse:
                    read.flag += -16
                else:
                    read.flag += 16
            outbam.write(read)

    inbam.close()
    outbam.close()


if __name__ == "__main__":
    main()