import collections
import sys
import click

"""
wirte by Jia Jinbu 2020.09.12

See more by --help.
"""

@click.command()
@click.option('-i', '--inreadinfo', help=('Input Read info file generated by merge_read_info.R or extract_read_info.py'
                                          'Required columns: "read_id", "chr", "read_start", "read_end", "mRNA", "mRNA_intron_num", '
                                          '"rel_gene_pos3", "f_feature_num", "l_feature_num", "retention_introns"'), 
                    required=True, type=click.Path(exists=True))
@click.option('-b', '--inbed', help="""Input bed file. The positions of exon and intron of represent transcripts.
                                       For example:


                        \b
                         chr1    3630    3913    AT1G01010.1_exon1       .       +
                         chr1    3913    3995    AT1G01010.1_intron1     .       +
                         chr1    3995    4276    AT1G01010.1_exon2       .       +
                         chr1    4276    4485    AT1G01010.1_intron2     .       +
                         chr1    4485    4605    AT1G01010.1_exon3       .       +
 
                         Note: the script will use "AT1G01010.1_exon1".split("_")[0] to extract the mRNA 
                         name, and use "AT1G01010.1_exon1".split("_")[1].startwith("e) check whether this 
                         feature is exon or intron. The number 1, 2, is the ordered from 5' mRNA to 3' mRNA.
                         """, 
                    required=True, type=click.Path(exists=True))                    
@click.option('-o', '--out', help="""Output Intron rel pos. Tab-seperated.                    

                        \b
                        intron_id           total_intron_num intron_num  rel_pos retention               read_core_id
                        AT1G01020.1_intron1     8                1        1439       0       7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687
                        AT1G01020.1_intron2     8                2        1300       0       7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687

                        rel_pos: the position of each intron relative to the 3' end of mRNA. 
                            (0 inidcate is just at the 3' end. thus the poistion of intron end is 
                            100, the 3' end position of mRNA is 110, the value is -10)

                        \b
                        retention: 0/1  0: spliced. 1: unspliced
                        """, 
                    required=True)
def main(inreadinfo, inbed, out):
    intron_pos3_data = read_intron_pos3(inbed)
    with open(out, 'w') as o:
        output_header = "intron_id\ttotal_intron_num\tintron_num\trel_pos\tretention\tread_core_id\n"
        o.write(output_header)
        select_columns = ["read_core_id", "chr", "read_start", "read_end", "mRNA", "mRNA_intron_num", "rel_mRNA_pos3", "f_feature_num", "l_feature_num", "retention_introns"]
        for read_id, chr_, read_start, read_end,  mRNA, mRNA_intron_num, rel_mRNA_pos3, \
                f_feature_num, l_feature_num, retention_introns in read_by_column_names(inreadinfo, select_columns):
            mRNA_intron_num, rel_gene_pos3, f_feature_num, l_feature_num = int(mRNA_intron_num), int(rel_mRNA_pos3), int(f_feature_num), int(l_feature_num)
            #1. remove mRNAs without intron
            if mRNA_intron_num == 0: continue
            #2. convert retention_introns to set. "3:4:5" to set([3, 4, 5]) 
            if not retention_introns or (retention_introns=="NA"):
                retention_introns = set()
            else:
                retention_introns = set([int(i) for i in retention_introns.split(":")])
            intron_datas = intron_pos3_data[mRNA]
            #3. iter each introns having been completely transcribed
            #The first intron is f_feature_num, and the last intron is l_feature_num-1, 
            #whatever the f_feature_type or l_feature_type is exon or intron.
            for intron_num in range(f_feature_num, l_feature_num):
                intron_pos3, intron_id = intron_datas[intron_num-1]
                retention = 1 if intron_num in retention_introns else 0
                rel_pos = rel_gene_pos3 - intron_pos3
                r = [intron_id, mRNA_intron_num, intron_num, rel_pos, retention, read_id]
                o.write("\t".join([str(s) for s in r]) + "\n")

def read_by_column_names(filein, select_column_names, sep="\t"):
    
    """
    Input:
    filein: tab-seperated file. The first line is column names.
    select_column_names: the list of names of column which you want to select.
    
    Output:
    A iterator.
    iter read each line of filein, and extract the values of select columns:
    
    For example:
    for column1, column2 in read_by_column_names(filein, ["column1", "column2"]):
        pass
    """
    
    def get_header(filein):
        header = ""
        with open(filein) as f:
            header = next(f)
        return header

    def get_header_index(filein, select_column_names, sep="\t"):
        select_indexs = []
        column_names = get_header(filein).rstrip("\n").split(sep)
        select_indexs = [column_names.index(column_name) for column_name in select_column_names]
        return select_indexs
    
    select_indexs = get_header_index(filein, select_column_names, sep)
    with open(filein) as f:
        next(f)
        for l in f:
            d = l.rstrip("\n").split(sep)
            yield([d[i] for i in select_indexs]) 

def read_intron_pos3(file_intron_exon):
    
    """
    Read exon intron bed file, extract the position of each intron relative to the 
    3' end of mRNA. (0 inidcate is just at the 3' end. thus the poistion of intron end is 
    100, the 3' end position of mRNA is 110, the value is -10)
    
    Input:
    file_exon_intron exon_intron bed file (0-based), 
            chr1    3630    3913    AT1G01010.1_exon1       .       +
            chr1    3913    3995    AT1G01010.1_intron1     .       +
    
    Output:
    A dict.
    Key: mRNA name  (AT1G01010.1)
    Value: a list of intron information, each element is [pos3, intron_id]. 
           ordered from the first intron to the last intron from mRNA 5' to mRNA 3'.
           intron_id: AT1G01010.1_intron1
           pos3: the position of intron relative to the 3' end of mRNA.
    """
    
    all_intron_poss = []
    mRNA2exons = collections.defaultdict(list) 
    for l in open(file_intron_exon):
        #AT1G01010.1_exon1
        #AT1G01010.1_intron1
        d = l.rstrip("\n").split("\t")
        exon_or_intron_name = d[3]
        chr_name = d[0]
        start, end = int(d[1])+1, int(d[2])
        strand = d[5]
        mRNA, exon_or_intron_id = exon_or_intron_name.split("_")
        if exon_or_intron_id.startswith("i"):
            all_intron_poss.append((exon_or_intron_name, mRNA, start, end))
        else:
            mRNA2exons[mRNA].append((chr_name, start, end, strand))
    mRNA_pos = {}
    for mRNA, mRNA_exons in mRNA2exons.items():
        total_exon_num = len(mRNA_exons)
        chr_name = mRNA_exons[0][0]
        strand = mRNA_exons[0][3]
        mRNA_start = min([i[1] for i in mRNA_exons])
        mRNA_end = max([i[2] for i in mRNA_exons])
        mRNA_pos[mRNA] = [chr_name, mRNA_start, mRNA_end, strand]
    intron_rel_poss = collections.defaultdict(list)
    for intron_id, mRNA, intron_start, intron_end in all_intron_poss:
        chr_name, mRNA_start, mRNA_end, strand = mRNA_pos[mRNA]
        if strand == "+":
            pos3 = int(intron_end) - mRNA_end
        else:
            pos3 = mRNA_start - int(intron_start)
        intron_rel_poss[mRNA].append([pos3, intron_id])
    for mRNA, d in intron_rel_poss.items():
        d.sort(key=lambda x: x[0])
    return intron_rel_poss
    
if __name__ == "__main__":
    main()


    
    