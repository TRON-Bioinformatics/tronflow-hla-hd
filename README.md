![TRON logo](https://tron-mainz.de/wp-content/uploads/2020/07/TRON_Logo_Science.svg "TRON logo")

--------

# TronFlow HLA typing pipeline

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Powered by Nextflow](https://img.shields.io/badge/powered%20by-Nextflow-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.nextflow.io/)



Nextflow (Di Tommaso, 2017) pipeline for HLA typing using HLA-HD (Kawaguchi, 2017).


## How to run it

Prepare an input table with the FASTQs for each sample with three tab-separated columns **without a header**

| Sample name          | FASTQ1                      | FASTQ2                  |
|----------------------|---------------------------------|------------------------------|
| sample_1             | /path/to/sample_1.1.fq      |    /path/to/sample_1.2.fq   |
| sample_2             | /path/to/sample_2.1.fq      |    /path/to/sample_2.2.fq  |


Run as indicated below.

```
$ nextflow run tron-bioinformatics/tronflow-hla-hd --help

N E X T F L O W  ~  version 19.07.0
Launching `main.nf` [intergalactic_shannon] - revision: e707c77d7b

Usage:
    nextflow run main.nf --input_files input_files --output output_folder

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name, FASTQ 1 and FASTQ 2
    The input file does not have header!
    Example input file:
    name1       fastq1.fq.gz    fastq2.fq.gz
    name2       fastq1.fq.gz    fastq2.fq.gz
    * output: output folder where results will be stored

Optional input:
    * read_length: the read length (default: 50)
    * hlahd_folder: the HLA-HD folder (default: /code/hlahd.1.2.0.1)
    * bowtie2_folder: the bowtie2 folder (default: /code/bowtie/2.3.4.3)
    * bowtie2_module: the module to load with bowtie2
    * ld_library_path: the value to set in LD_LIBRARY_PATH
    * cpus: the number of CPUs per sample (default: 15)
    * memory: the amount of memory per sample (default: 30g)

```


## References

* Kawaguchi S, Higasa K, Shimizu M, Yamada R, Matsuda F. HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. Hum Mutat. 2017 Jul;38(7):788-797. doi: 10.1002/humu.23230 Add to Citavi project by DOI. Epub 2017 May 12. PMID: 28419628
* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. 10.1038/nbt.3820
