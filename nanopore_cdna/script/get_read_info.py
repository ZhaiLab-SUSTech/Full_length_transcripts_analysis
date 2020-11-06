#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-20 14:27:41
@LastEditTime : 2020-05-29 18:47:46
@Description  : 将read比对到基因序列的信息整合到read_info中
                相当于注释reads
'''


import pickle
import pandas as pd
import numpy as np
from collections import defaultdict
import click

import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    )


def read_info_to_dict(df):
    logging.info('Start coverting read info to dict')

    df.set_index(['read_id'], inplace=True)
    df = df.loc[:, ['3\'adapter_start', 'tail_length', 'gene_id']].dropna()  # 未被注释的read直接drop
    df['gene_id'] = df['gene_id'].map(lambda x: x.split(','))
    df_dict = df.to_dict(orient='index')

    logging.info('Coverting read info to dict done!')
    return df_dict


@click.command()
@click.option('--paf_file', required=True)
@click.option('--polya_file', required=True)
@click.option('--outfile', required=True)
def main(paf_file, polya_file, outfile):

    '''
    1️⃣ 读取paf文件
    '''
    logging.info('Load paf file start')
    MM2_PAF_NAME = ['qname', 'qlen', 'qstart', 'qend', 
                'strand', 'tname', 'tlen', 
                'tstart', 'tend', 'residue_matches', 
                'alignment_block_length', 'mapq']
    mm2_paf_df = pd.read_csv(paf_file, sep='\t', names=MM2_PAF_NAME, header=None, usecols=[i for i in range(12)])
    logging.info('Load paf file end')

    '''
    2️⃣ 获取read_id对应的gene_id, 储存为字典
    '''
    logging.info('Get gene_id start')
    mm2_result = defaultdict(lambda: set())
    current_qname = None
    for item in mm2_paf_df.itertuples():
        if current_qname is None:
            current_qname = item.qname
        if current_qname != item.qname:
            mm2_result[current_qname] = ','.join(mm2_result[current_qname])
            current_qname = item.qname
        mm2_result[item.qname].add(item.tname)
    mm2_result[item.qname] = ','.join(mm2_result[current_qname])
    logging.info('Get gene_id end')

    '''
    3️⃣ 获取polyA信息, 并将比对结果整合到其中，导出read_info
    '''
    logging.info('Get polyA info start')
    def get_gene_id(read_id):
        if read_id in mm2_result:
            return mm2_result[read_id]
            
    read_info = pd.read_csv(polya_file, sep='\t')
    read_info['gene_id'] = read_info['read_id'].map(get_gene_id)
    read_info.to_csv(outfile, sep='\t', index=False)

    # 将结果保存为pickle格式，减少后面程序debug时运算耗时
    read_info_dict = read_info_to_dict(read_info)
    with open(outfile+'.pkl', 'wb') as o:
        pickle.dump(read_info_dict, o)
        
    logging.info('Get polyA info end')


if __name__ == "__main__":
    main()