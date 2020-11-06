'''
A pipline for nanopore cdna
You can run like this:
snakemake -j 6 -c "qsub -N {rulename} -l nodes=1:ppn={threads} -l mem=32g -j oe -l walltime=900:00:00"
'''

import os
import glob

configfile: 'config.yml'


# create log dir
if not os.path.exists('log'):
    os.mkdir('log')


SAMPLE_NAME = []
with open('basecalled_data/sample_name.txt', 'r') as f:
    for line in f:
        name = line.rstrip()
        SAMPLE_NAME.append(name)


rule all:
    input:
        expand('aligned_data/{sample_name}.adjust.sorted.tagged.bam', sample_name=SAMPLE_NAME),
        expand('polyadenylated_data/{sample_name}.full_len.bed.ir_stats.tsv', sample_name=SAMPLE_NAME),
        expand('polya_sites/{sample_name}.polya_cluster.bed', sample_name=SAMPLE_NAME)



# Adjust read strandedness and 3'adapter filter
rule adjust_strandedness:
    input:
        'basecalled_data/{sample_name}.fastq.gz',
    output:
        'basecalled_data/{sample_name}.adjust.fastq.gz'
    threads: 1
    params:
        polya_file=config['polya_data']
    shell:
        '''
python script/adjust_strandedness.py -i {input} -o {output} -p {params}
        '''


rule minimap2_map_genome:
    input:
        reads='basecalled_data/{sample_name}.adjust.fastq.gz',
    output:
        bam='aligned_data/{sample_name}.adjust.sorted.bam'
    threads: 64
    params:
        ref=config['genome'],
        bed=config['ann_bed'],
    shell:
        '''
minimap2 -t {threads} -ax splice -G 10000 --secondary=no --cs \
--junc-bed {params.bed} --sam-hit-only {params.ref} {input.reads} \
| samtools sort -@ {threads} -O bam -o {output.bam} -
samtools index -@10 {output}
        '''
        

rule bedtools_intersect:
    input:
        'aligned_data/{sample_name}.adjust.sorted.bam'
    output:
        temp('aligned_data/{sample_name}.adjust.sorted.bed.intersect')
    threads: 1
    params:
        config['ann_bed']
    shell:
        '''
bedtools intersect -abam {input} -b {params} -wo -s -split -bed > {output}
        '''


rule parse_bedtools_output:
    input:
        'aligned_data/{sample_name}.adjust.sorted.bed.intersect'
    output:
        'aligned_data/{sample_name}.adjust.sorted.bed.pkl'
    threads: 1
    shell:
        '''
python script/parse_bedtools_output.py -i {input} -o {output}
        '''


# add polya length and gene_id to bam tags
rule add_tags_to_bam:
    input:
        infile='aligned_data/{sample_name}.adjust.sorted.bam',
        bed_intersect='aligned_data/{sample_name}.adjust.sorted.bed.pkl'
    output:
        'aligned_data/{sample_name}.adjust.sorted.tagged.bam'
    params:
        config['polya_data']
    threads: 1
    shell:
        '''
python script/add_tag_to_bam.py -i {input.infile} -o {output} --read_info {params} --bed_intersect {input.bed_intersect}
samtools index -@10 {output}
        '''


rule get_polyadenylated_reads:
    input:
        'aligned_data/{sample_name}.adjust.sorted.tagged.bam'
    output:
        'polyadenylated_data/{sample_name}.polyadenylated.bam'
    threads: 1
    shell:
        '''
python script/get_polyadenylated_reads.py -i {input} -o {output}
samtools index -@20 {output}
        '''


rule get_polya_cluster:
    input:
        bam='polyadenylated_data/{sample_name}.polyadenylated.bam'
    output:
        'polya_sites/{sample_name}.polya_cluster.bed',
        'polya_sites/{sample_name}.polya_cluster.summit.bed',
        'polya_sites/{sample_name}.major_polya_cluster.bed',
        'polya_sites/{sample_name}.major_polya_cluster_summit.bed'
    params:
        suffix='polya_sites/{sample_name}',
        gene_bed='/public/home/mowp/db/Arabidopsis_thaliana/isoform/araport11.gene.bed'
    threads: 32
    shell:
        '''
python script/get_polya_cluster.py --infile {input.bam} --gene_bed {params.gene_bed} --out_suffix {params.suffix} -t {threads}
        '''

rule get_full_len_transcripts:
    input:
        'polyadenylated_data/{sample_name}.polyadenylated.bam'
    output:
        'polyadenylated_data/{sample_name}.full_len.bam'
    threads: 1
    params:
        first_exon_path='/public/home/mowp/db/Arabidopsis_thaliana/representative_gene_model/representative_gene_first_exon.bed'
    shell:
        '''
python script/remove_five_prime_lost.py -i {input} -o {output} --first_exon_path {params.first_exon_path}
sleep 15
samtools index -@ 16 {output}
        '''


rule bedtools_intersect_for_splice_stats:
    input:
        'polyadenylated_data/{sample_name}.full_len.bam'
    output:
        'polyadenylated_data/{sample_name}.full_len.bed.intersect'
    params:
        repr_bed='~/db/Arabidopsis_thaliana/representative_gene_model/araport11.representative.gene_model.bed'
    shell:
        '''
bedtools intersect -abam {input} -b {params.repr_bed} -wo -s -split -bed > {output}
        '''


rule get_splice_stats:
    input:
        'polyadenylated_data/{sample_name}.full_len.bed.intersect'
    output:
        ir='polyadenylated_data/{sample_name}.full_len.bed.ir_stats.tsv',
        es='polyadenylated_data/{sample_name}.full_len.bed.es_stats.tsv'
    shell:
        '''
python script/get_splice_stats.py -i {input} --ir_out {output.ir} --es_out {output.es}
        '''


STRAND = {'plus': '+', 'minus': '-'}
rule bamtobigwig_split:
    input:
        '{dirname}/{sample_name}.bam'
    output:
        bedgraph=temp('{dirname}/{sample_name}.cov.{strand}.bedgraph'),
        bw='{dirname}/{sample_name}.cov.{strand}.bw'
    threads: 1
    params:
        strand='{strand}',
        chrom_size='/public/home/mowp/db/Arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.chrom.sizes'
    shell:
        '''
bedtools genomecov -split -strand {params.strand} -bg -ibam {input} > {output.bedgraph}
bedGraphToBigWig {output.bedgraph} {params.chrom_size} {output.bw}
        '''


rule bamtobigwig_3end:
    input:
        '{dirname}/{sample_name}.bam'
    output:
        bedgraph=temp('{dirname}/{sample_name}.3end.{strand}.bedgraph'),
        bw='{dirname}/{sample_name}.3end.{strand}.bw'
    threads: 1
    params:
        strand='{strand}',
        chrom_size='/public/home/mowp/db/Arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.chrom.sizes'
    shell:
        '''
bedtools genomecov -3 -strand {params.strand} -bg -ibam {input} > {output.bedgraph}
bedGraphToBigWig {output.bedgraph} {params.chrom_size} {output.bw}
        '''
