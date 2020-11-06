#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-24 19:07:01
@LastEditTime : 2020-08-01 15:06:23
@Description  : 
'''


from collections import Counter
import numpy as np
import pysam
import click
import pyranges as pr
from concurrent.futures import ProcessPoolExecutor


def get_entropy(site_counter, total_count):
    '''
    计算这个PAC内的polya site的entropy
    '''
    entropy = 0
    for site in site_counter:
        p = site_counter[site]/total_count
        entropy += -p*np.log2(p)
    return entropy



STRAND_TO_BOOL = {'-': True, '+': False}


def get_three_end(infile, gene_id, gene_model):
    '''
    获取polyadenylated read 3’end 位置
    '''
    chrom, start, end, _, strand = gene_model.loc[gene_id, :].values
    strand = STRAND_TO_BOOL[strand]
    read_count = 0
    polya_sites = []
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        for read in inbam.fetch(chrom, start, end):
            # 取该位置上与基因方向一致的reads
            if strand is read.is_reverse:
                polya_len = read.get_tag('pa')
                read_gene_id = read.get_tag('gi')
                if polya_len >= 15 and read_gene_id in {'None', gene_id}:
                    read_count += 1
                    if not read.is_reverse:
                        polya_sites.append(read.reference_end)
                    else:
                        polya_sites.append(read.reference_start*-1)
    
    # 选择表达量高的基因
    if read_count < 10:
        return None
    
    polya_sites.sort()
    
    # 将24nt以内的点合并
    total_site_count = len(polya_sites)
    pac_list = None
    for polya_site in polya_sites:
        if pac_list is None:
            pac_list = [[polya_site]]
        else:
            if polya_site-pac_list[-1][-1] <= 24:
                pac_list[-1].append(polya_site)
            else:
                pac_list.append([polya_site])

    # 去除read数少于3的pac
    # 去除read比例少于基因区总read数1%的pac
    major_cluster_site_count = 0
    summit = []

    polya_cluster = []
    polya_cluster_summit = []

    polya_cluster_major = ''
    polya_cluster_summit_major = ''

    polya_cluster_last = ''
    polya_cluster_summit_last = ''

    *_, strand = gene_model.loc[gene_id, :].values
    n = 0
    for pac in pac_list:
        # 去除read数少于3的pac
        # 去除read比例少于基因区总read数1%的pac
        if len(pac) < 3 or len(pac)/total_site_count < .1:
            continue

        site_counter = Counter(pac)
        site_counter_most = site_counter.most_common(3)

        if site_counter_most[0][1] >= 3:
            n += 1
            start = min(abs(pac[0]), abs(pac[-1]))
            end = max(abs(pac[0]), abs(pac[-1]))
            entropy = get_entropy(site_counter, len(pac))
            polya_cluster.append(f'{chrom}\t{start-1}\t{end}\t{gene_id}_{n}\t{len(pac)}\t{strand}\t{entropy:.3f}\n')

            # 计算cluster的峰值
            summit = abs(site_counter_most[0][0])
            # 最后一列为summit与pac中所有counts的比值
            # 比值越高，summit越显著
            polya_cluster_summit.append(f'{chrom}\t{summit-1}\t{summit}\t{gene_id}_{n}\t{site_counter_most[0][1]}\t{strand}\t{site_counter_most[0][1]/len(pac):.3f}\n')

            if major_cluster_site_count < len(pac):
                major_cluster_site_count = len(pac)
                polya_cluster_major = polya_cluster[-1]
                polya_cluster_summit_major = polya_cluster_summit[-1]
    
    if len(polya_cluster) == 0:
        return None
    
    polya_cluster_last = polya_cluster[-1]
    polya_cluster_summit_last = polya_cluster_summit[-1]
    
    return polya_cluster, polya_cluster_summit, polya_cluster_major, polya_cluster_summit_major, polya_cluster_last, polya_cluster_summit_last


@click.command()
@click.option('--infile', required=True)
@click.option('--gene_bed', required=True)
@click.option('--out_suffix', required=True)
@click.option('-t', '--threads', required=False, default=10)
def main(infile, gene_bed, out_suffix, threads):
    # 读取gene_bed信息
    gene_model = pr.read_bed(gene_bed, as_df=True)
    gene_model = gene_model.set_index(['Name'])

    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        for gene_id in gene_model.index:
            if gene_model.at[gene_id, 'Chromosome'] not in {'Mt', 'Pt'}:
                results.append(e.submit(get_three_end, infile, gene_id, gene_model))
    

    o_polya_cluster = open(f'{out_suffix}.polya_cluster.bed', 'w')
    o_polya_cluster_summit = open(f'{out_suffix}.polya_cluster.summit.bed', 'w')
    o_major_polya_cluster = open(f'{out_suffix}.major_polya_cluster.bed', 'w')
    o_major_polya_cluster_summit = open(f'{out_suffix}.major_polya_cluster_summit.bed', 'w')
    o_last_polya_cluster = open(f'{out_suffix}.last_polya_cluster.bed', 'w')
    o_last_polya_cluster_summit = open(f'{out_suffix}.last_polya_cluster_summit.bed', 'w')

    for res in results:
        result = res.result()
        if result is not None:
            polya_cluster, polya_cluster_summit, polya_cluster_major, polya_cluster_summit_major, polya_cluster_last, polya_cluster_summit_last = result
            for item in polya_cluster:
                o_polya_cluster.write(item)
            for item in polya_cluster_summit:
                o_polya_cluster_summit.write(item)
                
            o_major_polya_cluster.write(polya_cluster_major)
            o_major_polya_cluster_summit.write(polya_cluster_summit_major)
            o_last_polya_cluster.write(polya_cluster_last)
            o_last_polya_cluster_summit.write(polya_cluster_summit_last)
            
    o_polya_cluster.close()
    o_polya_cluster_summit.close()
    o_major_polya_cluster.close()
    o_major_polya_cluster_summit.close()
    o_last_polya_cluster.close()
    o_last_polya_cluster_summit.close()


if __name__ == "__main__":
    main()
    