#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-21 11:00:38
@LastEditTime : 2020-05-21 11:30:59
@Description  : 
'''


import click
import pysam


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
def main(infile, outfile):

    with pysam.AlignmentFile(infile, 'rb') as inbam, \
        pysam.AlignmentFile(outfile, 'wb', template=inbam) as outbam:
        for read in inbam:
            '''
            polyadenylated reads:
                1. polyA tail length >= 15
            '''
            if read.get_tag('pa') >= 15:
                outbam.write(read)


if __name__ == "__main__":
    main()