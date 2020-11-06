#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-20 19:22:57
@LastEditTime : 2020-06-24 16:51:09
@Description  : get elongating reads
'''


import pysam
import pandas as pd
import numpy as np
import click


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
def main(infile, outfile):
    with pysam.AlignmentFile(infile, 'rb') as inbam, \
        pysam.AlignmentFile(outfile, 'wb', template=inbam) as outbam:
        for read in inbam:
            '''
            elongating reads:
                1. polyA tail length < 15
                2. distance between the start position of the 3'adapter 
                and the end position of the genome mapping region <= 5
            '''
            if read.get_tag('pa') < 15 and read.get_tag('ga') <= 5:
                outbam.write(read)


if __name__ == "__main__":
    main()