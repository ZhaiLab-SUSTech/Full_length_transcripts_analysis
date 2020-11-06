import pysam
import sys

"""
python this.py file_bam fileout
根据转录组比对文件，提取每个read的exon信息。由于read来源于RNA，有可能有剪切，根据剪切情况，区分不同的exon。
这里忽略缺失和删除，只用skip信息。返回的exon起始和终止位置是1-based的，而且终止位置包含在exon内。
由于一个read有可能比对到多个位置，因此这里将read名添加上了位置比对信息。

read名重命名为（逗号隔开的字符串）:
read_origin_name,chr_name,start,end,strand,total_exon_num
start和end是比对位置，1-based，end包含在read内。
total_exon_num是该read总共有多少个exon。

输出文件格式(tab隔开):
每一行是一个exon，名字为read_name加exon_num(按染色体位置先后排列，不考虑比对方向)（逗号隔开）。
read_name,exon_num chr_name start end strand
"""

file_bam, fileout = sys.argv[1:3]

SEP = ","

#file_bam = "all.insert.minimap2.bam"
#fileout = "all.insert.minimap2.read_exon.txt"
bam_obj = pysam.AlignmentFile(file_bam, "rb")

def get_read_blocks(read):
    '''
    这里忽略I和D，根据N（skipped region）来划分exon。
    返回exon的列表，每个元素是一个exon的起始和终止位置，注意都是1-based。
    而且终止位置包含在exon内。
    
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    '''
    blocks = []
    start = read.reference_start + 1
    end = start - 1
    for (flag, length) in read.cigartuples:
        if flag == 4 or flag == 5 or flag == 1: continue
        if flag == 0 or flag == 2:
            end += length
        if flag == 3:
            blocks.append([start, end])
            start = end + length + 1
            end = start - 1
    blocks.append([start, end])
    return blocks

def get_introns(blocks):
    introns = []
    if len(blocks) > 1:
        last_start = blocks[0][1] + 1
        for start, end in blocks[1:]:
            introns.append([last_start, start-1])
            last_start = end + 1
    return introns

def get_read_introns(read):
    return get_introns(get_read_blocks(read))

with open(fileout, 'w') as o:
    for read in bam_obj.fetch():
        #reference_id, reference_end, reference_start, 
        read_strand = "-" if read.is_reverse else "+"        
        read_blocks = get_read_blocks(read)
        total_block_num = len(read_blocks)
        read_name = SEP.join([read.query_name,
                              read.reference_name, 
                              str(read.reference_start + 1), 
                              str(read.reference_end),
                              read_strand,
                              str(total_block_num)])
        for i, (start, end) in enumerate(read_blocks):
            #exon_num = i + 1 if read_strand == "+" else total_block_num - i
            exon_num = i + 1
            o.write("\t".join([read_name + SEP + str(exon_num), read.reference_name, str(start), str(end), read_strand]) + "\n")

    
    