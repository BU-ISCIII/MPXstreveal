# MPXstreveal

### Overview

 &emsp;  &emsp; Monkeypox STR reveal (MPXstreveal) a new draft approach for detecting and curating low complexity regions in monkeypox genomes using massive parallel sequencing. STR search algorithm is based in the [STRsearch pipeline](https://github.com/AnJingwd/STRsearch) [[1]](https://hereditasjournal.biomedcentral.com/articles/10.1186/s41065-020-00120-6)

<!---
 TODO ### Algorithm description: 

 &emsp;  &emsp; Briefly, STRsearch employs an iterative algorithm to obtain the longest continuous interval composed by all motifs of STR sequence structure without a priori assumptions on allele size. The actual STR region is determined by comparing the position of repeat patterns with the best matching location of flanking sequences in reads. Ultimately, allele size is calculated not only for repeat patterns, but also indels that are actually in the STR region.  


TODO ### Installation

To obtain STRsearch, use:

```Git
git clone https://github.com/AnJingwd/STRsearch.git
```

or

```shell
wget https://github.com/AnJingwd/STRsearch/archive/master.zip
```

TODO ### Prerequisite

 &emsp;  &emsp; The following linux utilities are needed and the full path of  them on your  local machine  should be provided in conf.py file

1. bwa  (v1.7 or higher) 
2. samtools  (v1.7 or higher) 
3. bamToFastq (v2.17.0 or higher)
4. seqtk  (v1.2 or higher)
5. usearch (v11 or higher)  [Download](https://www.drive5.com/usearch/download.html)

Additionally, the following Python modules are required.

1. numpy
2. argparse
3. pathlib



TODO ### Configuration file format:

 &emsp;  &emsp; The first step for  STR analysis with STRsearch is to create a configuration file with your custom set of STR loci. One way to do this is by referring to the most up-to-date revised forensic STR sequence guide  and a worksheet can be downloaded from [link](https://strider.online/bundles/strbaseclient/downloads/Forensic_STR_Sequence_Structure_Guide_v5.xlsx). You will need to make a configuration file with the following columns present: 

| Chr   | Start     | End       | Period | Reference allele | Marker  | STR       | STR sequence structure                                       | Stand | 5' Flanking  sequence | 3' Flanking  sequence |
| ----- | --------- | --------- | ------ | ---------------- | ------- | --------- | ------------------------------------------------------------ | ----- | --------------------- | --------------------- |
| chr1  | 7442891   | 7442934   | 4      | 11               | Marker1 | D1GATA113 | [GATA]n                                                      | +     | ACTTGCTTCCTAGAT       | TTCCTATAGCCTCAA       |
| chr21 | 20554291  | 20554417  | 4      | 29               | Marker2 | D21S11    | [TCTA]n  [TCTG]n [TCTA]n ta [TCTA]n tca [TCTA]n tccata [TCTA]n TA [TCTA]n | +     | CCAAGTGAATTGCCT       | TCGTCTATCTATCCA       |
| chrX  | 149710971 | 149711038 | 4      | 15               | Marker3 | DXS7423   | [TGGA]n aggacaga  [TGGA]n                                    | +     | AAATGAATGAGTATG       | TGGGGAGGAAATCTG       |
| chrY  | 15752608  | 15752715  | 3      | 27               | Marker4 | DYS612    | [CCT]n CTT [TCT]n  CCT [TCT]n                                | +     | AGGTTCAGAGGTTTG       | GTCACTTTTCCAAAT       |
| chrY  | 20842518  | 20842573  | 4      | 14               | Marker5 | DYS385a   | [TTTC]n                                                      | -     | TCCTTTCTTTTTCTC       | CCTTCCTTCCTTCCT       |

1. Column 1 : chromosome  (must)
2. Column 2 : start coordinate of the STR  (must)
3. Column 3  : end coordinate of the STR  (must)
4. Column 4  : period of the STR  (must)
5. Column 5 : reference copy number (option)
6. Column 6 : Marker name (option) 
7. Column 7 : STR name (option) 
8. Column 8 :   Reference Sequence repeat region sequence structure summary  (must)
9. Column 9 : stand ("+" means positive stand;"-" means negative stand)  (must)
10. Column 10 : 5' flanking sequence of repeat region  (must)
11. Column 11 : 3' flanking sequence of repeat region  (must)



Note some columns are not used. You can put any value in the non-required columns, just make sure there are at least 11 columns with the required information listed above. Importanly, flanking sequences are necessarily adjacent to STR repeat region.  



TODO ### Inputs

FASTQ file or BAM-file from singe-end or paird-end sequencing platforms

TODO ### Output

- genotypes.txt: genotypes on each targeted locus

- multiple_alleles.txt: all alleles identified on each targeted locus

- qc_matrix.txt:  a quality control matrix including several sequence properties (total bases, sequencing quality score, number of allocated reads, distance distribution of STR repeat sequence to end of reads, allele read depth)  


TODO ### Usage examples

**1. run with  default parameters**

for paired-end sequencing

```shell
python3 pipeline.py from_fastq \
--working_path example/test_results/ \
--sample test \
--fq1 example/test_data/test_R1.fastq \
--fq2 example/test_data/test_R2.fastq \
--ref ucsc.hg19.fasta
```



```shell
python3 pipeline.py from_bam \
--working_path example/test_results \
--sample test \
--sex male \
--bam example/test_results/alignments/test.bam \
--ref_bed example/ref_test.bed \
--genotypes example/test_results/test_genotypes.txt \
--multiple_alleles example/test_results/test_multiple_alleles.txt \
--qc_matrix example/test_results/test_qc_matrix.txt
```



for single-end sequencing

```shell
python3 pipeline.py  \
--type single \
from_fastq \
--working_path example/test_results/ \
--sample test \
--fq1 example/test_data/test_R1.fastq \
--ref ucsc.hg19.fasta
```



```shell
python3 pipeline.py \
--type single \
from_bam \
--working_path example/test_results \
--sample test \
--sex male \
--bam example/test_results/alignments/test.bam \
--ref_bed example/ref_test.bed \
--genotypes example/test_results/test_genotypes.txt \
--multiple_alleles example/test_results/test_multiple_alleles.txt \
--qc_matrix example/test_results/test_qc_matrix.txt
```



**2. run with self-defined parameters**

```shell
python3 pipeline.py \
--assemble_pairs True \
--reads_threshold 50 \
--stutter_ratio 0.6 \
--num_threads 8 \
--num_processors 8 \
from_bam \
--working_path example/test_results \
--sample test \
--sex male \
--bam example/test_results/alignments/test.bam \
--ref_bed example/ref_test.bed \
--genotypes example/test_results/test_genotypes.txt \
--multiple_alleles example/test_results/test_multiple_alleles.txt \
--qc_matrix example/test_results/test_qc_matrix.txt
```



TODO ### Options

**Default parameters**

| Option            | Value Type | Default | Summary                                          |
| ----------------- | ---------- | ------- | ------------------------------------------------ |
| --help            |            | false   | display the help message                         |
| --type            | str        | paired  | (option) The sequencing type                     |
| --assemble_pairs  | bool       | False   | (option) if True, paired-end reads are assembled |
| --reads_threshold | int        | 30      | (option) The analytical threshold for reads      |
| --stutter_ratio   | float      | 0.5     | (option) The stutter ratio                       |
| --num_threads     | int        | 4       | (option) The number of multiple threads          |
| --num_processors  | int        | 4       | (option) The number of multiprocess              |

**Sub command**

1. from_bam

| Option             | Value Type | Default | Summary                                      |
| ------------------ | ---------- | ------- | -------------------------------------------- |
| --help             |            | false   | display the help message                     |
| --working_path     | str        | null    | (must) The working path                      |
| --sample           | str        | null    | (must) The sample name                       |
| --sex              | str        | null    | (must) The sample sex                        |
| --bam              | str        | null    | (must) The input BAM-file                    |
| --ref_bed          | str        | null    | (must) The configuration file of STRs        |
| --genotypes        | str        | null    | (must) The output for STR genotypes          |
| --multiple_alleles | str        | null    | (must) The output for multiple alleles       |
| --qc_matrix        | str        | null    | (must) The output for quality control matrix |



2. from_fastq


| Option         | Value Type | Default | Summary                                                      |
| -------------- | ---------- | ------- | ------------------------------------------------------------ |
| --help         |            | false   | display the help message                                     |
| --working_path | str        | null    | (must) The working path                                      |
| --sample       | str        | null    | (must) The sample name                                       |
| --fq1          | str        | null    | (must) The in1.fq                                            |
| --fq2          | str        | null    | (option) The in2.fq for paired-end sequencing                |
| --ref          | str        | null    | (must) The reference genome fasta and index file in the same path |



TODO ### Run with Docker

To obtain STRsearch Docker image, use:

```GO
docker pull anjing123/strsearch:latest
```



```
LOCAL_PATH/app/
├── ref
│   ├── ucsc.hg19.fasta
│   ├── ucsc.hg19.fasta.amb
│   ├── ucsc.hg19.fasta.ann
│   ├── ucsc.hg19.fasta.bwt
│   ├── ucsc.hg19.fasta.pac
│   └── ucsc.hg19.fasta.sa
├── ref_test.bed
└── test_data
    ├── test.bam
    ├── test_R1.fastq
    └── test_R2.fastq
```



```shell
docker run -v LOCAL_PATH/app/:/app/ -w /app/ -it anjing123/strsearch:latest from_fastq \
--working_path /app/test_results/ \
--sample test \
--fq1 /app/test_data/test_R1.fastq \
--fq2 /app/test_data/test_R2.fastq \
--ref /app/ref/ucsc.hg19.fasta
```



```shell
docker run -v LOCAL_PATH/app/:/app/ -w /app/ -it anjing123/strsearch:latest from_bam \
--working_path /app/test_results/ \
--sample test \
--sex male \
--bam /app/test_results/alignments/test.bam \
--ref_bed /app/ref_test.bed \
--genotypes /app/test_results/test_genotypes.txt \
--multiple_alleles /app/test_results/test_multiple_alleles.txt \
--qc_matrix /app/test_results/test_qc_matrix.txt
```
-->

**This repo is a work in progress that woul soon be updated**

### Acknowledgements

Thanks to the team that developed [STRsearch pipeline](https://github.com/AnJingwd/STRsearch) in which is based our code with slight modifications.


### References

1. [Wang, D., Tao, R., Li, Z. et al. STRsearch: a new pipeline for targeted profiling of short tandem repeats in massively parallel sequencing data. Hereditas 157, 8 (2020). https://doi.org/10.1186/s41065-020-00120-6](https://hereditasjournal.biomedcentral.com/articles/10.1186/s41065-020-00120-6)

<!---
TODO ### Contact

Sara Monzón 
Sarai Varona
Isabel Cuesta
-->
