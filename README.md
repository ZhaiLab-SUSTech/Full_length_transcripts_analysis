
<!--
 * @Author       : windz
 * @Date         : 2020-11-10 21:49:47
 * @LastEditTime : 2020-11-10 22:11:40
 * @Description  : README
-->


# Full-length transcript analysis 

This repository is the pipeline and scripts of analysing the Nanopore and PacBio datasets for genome-wide characterizing isoform-specific poly(A) tail length, and includes snakemake workflow and ipython notebooks for generating the paper figures.

## Software Dependencies

The snakemake workflow is implemented in Python3 and requires [minimap2](https://github.com/lh3/minimap2) and [bedtools](https://bedtools.readthedocs.io/) as addtional software. Dependent Packages are listed below:

```
pysam
numpy
pandas
matplotlib
pyranges
click
scipy
```