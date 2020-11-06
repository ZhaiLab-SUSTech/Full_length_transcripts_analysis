#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-14 10:46:37
LastEditTime : 2020-10-26 11:06:30
@Description  : add polya length and remove wrong strandness reads
                tags:
                    'ga': distance between the start position of the 3’adapter 
                          and the end position of the genome mapping region
                    'pa': polya tail length
                    'gi': gene_id
'''


import concurrent.futures
import pickle
import pysam
import pandas as pd
import numpy as np
import click


import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    )


def load_read_info(read_info):
    '''
    load read info from polyaCaller result

    Return:
        {'3da1bb92-be52-47c6-affd-5ef02febc17b': 
            {"3'adapter_start": 69.0,
            'tail_length': 0.99},
        'dfdca161-7d81-4838-b7cb-6a96c47d6d2c': 
            {"3'adapter_start": 551.0,
            'tail_length': 0.76},
        ...
        }
    '''
    logging.info('Load read info')

    USE_COLS = ['read_id', 'read_type', 'mapped_start', 'mapped_end', '3\'adapter_start', 'tail_length']
    df = pd.read_csv(read_info, sep='\t', index_col=['read_id'], usecols=USE_COLS)

    df['gap_size'] = float('nan')
    df['polya'] = df["3'adapter_start"] - df['mapped_end']
    df['polyt'] = df['mapped_start'] - df["3'adapter_start"] - 1
    df.loc[df['read_type'] == 'polyA', 'gap_size'] = df['polya'] 
    df.loc[df['read_type'] == 'polyT', 'gap_size'] = df['polyt'] 
    read_info_dict = df.loc[:, ['tail_length', 'gap_size']].to_dict(orient='index')
        
    logging.info('Load read info done!')
    return read_info_dict


def get_intersect_info(bed_intersect):
    logging.info('Load bed intersect info')
    with open(bed_intersect, 'rb') as f:
        intersect_info_dict = pickle.load(f)
    logging.info('Load bed intersect info done!')
    return intersect_info_dict


def find_introns(r):
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    res = []
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    base_position = r.reference_start
    for op, nt in r.cigartuples:
        if op in match_or_deletion:
            base_position += nt
        elif op == BAM_CREF_SKIP:
            junc_start = base_position
            base_position += nt
            res.append((junc_start, base_position))
    return res


def find_exons(r):
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    res = []
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    exon_start = r.reference_start
    length = 0
    for op, nt in r.cigartuples:
        if op in match_or_deletion:
            length += nt
        elif op == BAM_CREF_SKIP:
            res.append((exon_start, exon_start+length))
            exon_start = res[-1][1]+nt
            length = 0
    res.append((exon_start, exon_start+length))
    return res


def is_exceed_extend(read, gene_len):
    '''
    判断intron是否过长
    如果该intron比所在基因还长
    '''
    introns = np.array(find_introns(read))
    if len(introns) > 0:
        intron_len = introns[:, 1] - introns[:, 0]
        return (intron_len > gene_len).any()
    else:
        return False


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
@click.option('--read_info', required=True)  # 提供polya长度和3'adapter信息
@click.option('--bed_intersect', required=True)  # parse_bedtools_output的pickle_file
def main(infile, outfile, read_info, bed_intersect):
    logging.info('Waiting for all subprocesses done...')
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:
        read_info_dict = e.submit(load_read_info, read_info)
        intersect_info_dict = e.submit(get_intersect_info, bed_intersect)
    '''    
    # polyaCaller信息
    read_info:
        polyacaller结果, 记录read 3'adapter, polya_len
        key为read_id

    # 基因注释的信息
    intersect_info_dict:
        {'669710e6-e0ef-45c2-a5ed-8c70bfe9facd': {
            'cov': 311,
            'is_splicing_intermediates': False,
            'gene_len': 7705,
            'gene_id': 'AT1G01040.2'},
        ...}
    '''
    read_info_dict = read_info_dict.result()
    intersect_info_dict = intersect_info_dict.result()
    logging.info('All subprocesses done.')

    logging.info('Start main function')
    inbam = pysam.AlignmentFile(infile, 'rb') 
    outbam = pysam.AlignmentFile(outfile, 'wb', template=inbam)

    for read in inbam:
        if not read.is_unmapped \
            and not read.is_supplementary \
                and read.query_name in read_info_dict:
            
            # 获取read的polyA长度
            read_tail_length = read_info_dict[read.query_name]['tail_length']

            # 判断是否为剪切中间体
            try:
                if intersect_info_dict[read.query_name]['is_splicing_intermediates'] \
                    and read_tail_length < 15:
                    continue
            except KeyError:
                pass

            # 判断intron是否由于比对错误导致过长的情况
            try:
                if is_exceed_extend(read, intersect_info_dict[read.query_name]['gene_len']):
                    continue
            except KeyError:
                logging.info(f'read {read.query_name} not in bed intersect file')

            # 获取gap长度
            gap_size = read_info_dict[read.query_name]['gap_size']
            
            if gap_size < -1:
                continue

            # 过滤adapter找错情况，理论上gap长度就是polya长度
            # 如果gap长度比polyA还要长，有可能是接头找错
            if read_tail_length> 15 and gap_size > read_tail_length :
                continue

            # 获取gene_id
            try:
                gene_id = intersect_info_dict[read.query_name]['gene_id']
            except KeyError:
                gene_id = 'None'

            # 3'adapter到read的距离
            read.set_tag('ga', gap_size)
            read.set_tag('pa', read_tail_length)
            read.set_tag('gi', gene_id)
            outbam.write(read)

    
    inbam.close()
    outbam.close()


if __name__ == "__main__":
    main()
    