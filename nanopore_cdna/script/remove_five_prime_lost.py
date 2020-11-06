#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-16 15:37:46
@LastEditTime : 2020-06-24 16:58:34
@Description  : remove reads lost 5' end
                the 5' end is located downstream of the annotated first exon, 
                likely due to incomplete reverse transcription
'''


import pysam
import pandas as pd
import click


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
@click.option('--first_exon_path', required=True)
def main(infile, outfile, first_exon_path):

    first_exon_path_df = pd.read_csv(
        first_exon_path, 
        sep='\t',
        index_col=['gene_id']
    )

    with pysam.AlignmentFile(infile, 'rb') as inbam, \
        pysam.AlignmentFile(outfile, 'wb', template=inbam) as outbam:
        for read in inbam:
            read_gene_id = read.get_tag('gi')
            if read_gene_id == 'None':
                continue
            item = first_exon_path_df.loc[read_gene_id]
            is_intact = False
            
            '''
            single_exon基因筛选条加：
                * 要5'end位置在基因起始位置下游30%范围上游
                * 且read比对到基因组的长度（头减尾，非序列长度，
                  小于基因长度的1.5倍太长的可能就不是那个基因了
            '''
            if item['is_single_exon']:
                exon_len = item['exon_end']-item['exon_start']

                if (read.reference_end-read.reference_start)/exon_len >= 1.5:
                    continue
                
                if not read.is_reverse \
                    and read.reference_start < item['exon_start']+exon_len*.3:
                    is_intact = True
                elif read.is_reverse \
                    and read.reference_end > item['exon_end']-exon_len*.3:
                    is_intact = True
    
            elif not read.is_reverse:
                read_fivePrime_pos = read.reference_start
                first_exon_end = item['exon_end']
                if read_fivePrime_pos <= first_exon_end:
                    is_intact = True
                    
            else:
                read_threePrime_pos = read.reference_end
                first_exon_start = item['exon_start']
                if read_threePrime_pos >= first_exon_start:
                    is_intact = True
            
            if is_intact:
                outbam.write(read)


if __name__ == "__main__":
    main()