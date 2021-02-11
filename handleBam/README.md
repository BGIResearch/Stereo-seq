# HandleBam

## Introduction

HandleBam implement these steps in RNA-seq pipeline:
1. mapping quality filter
2. deduplication
3. set annotations
4. stat gene expression data

optional functions:

1. umi correction
2. sequencing saturation

## Compile

### Platform & Environment

* centos-7.0+
* gcc-9.1.0
* cmake-3.17.2

### Prerequisites

| Library    | Version | Description                        | Link                                                                        |
| ---------- | ------- | ---------------------------------- | --------------------------------------------------------------------------- |
| htslib     | 1.9.0   | process bam/sam data               | https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 |
| spdlog     | 1.5.0   | logging module                     | https://github.com/gabime/spdlog/archive/v1.5.0.zip                         |
| CLI11      | 1.9.0   | parse command line parameters      | https://github.com/CLIUtils/CLI11/releases/download/v1.9.0/CLI11.hpp        |
| libdeflate | 1.5     | accelerate bgzf IO                 | https://github.com/ebiggers/libdeflate/archive/v1.5.zip                     |
| ygg        | master  | interval tree for find overlapping | https://github.com/tinloaf/ygg.git                                          |
| doctest    | 2.3.7   | optional for unittest              | https://github.com/onqtam/doctest/archive/2.3.7.tar.gz                      |
| fftw       | 3.3.8   | for KDE                            |                                                                             |

### Compile

1. enter the code path
2. modify the value of `libPath` in file *script/build.sh*
3. run *script/build.sh*
4. the compile results are saved in this directory: *install/bin install/lib*

### UnitTest

Unit testing is already supported.  
Repeat steps in *Compile* and modify step three: `sh ./script/build.sh ON` , 
then the binary executable file named *unittest* is saved in *install/bin*

## Run

### Usage

```
HandleBam: mapping quality filter, deduplication, set annotation, stat gene expression data.
Usage: ./install/bin/handleBam [OPTIONS]

Options:
  -h,--help                             Print this help message and exit
  -I,-i TEXT REQUIRED                   Input bam filename or file list separated by comma
  -O,-o TEXT REQUIRED                   Output bam filename
  -A,-a TEXT:FILE REQUIRED              Input annotation filename
  -S,-s TEXT REQUIRED                   Output summary filename
  -E,-e TEXT REQUIRED                   Output barcode gene expression filename
  -Q,-q INT:POSITIVE                    Set mapping quality threshold, default 10
  -C,-c INT:POSITIVE                    Set cpu cores, default detect
  --save_lq                             Save low quality reads, default false
  --save_dup                            Save duplicate reads, default false
  --anno_mode INT:INT in [0 - 2]        Select annotation mode, default 2
                                        0->'v1'
                                        1->'v2'
                                        2->'v3'
  --umi_on                              Open umi correction, default off
  --umi_min_num INT:POSITIVE            Minimum umi number for correction, default 5
  --umi_mismatch INT:POSITIVE           Maximum mismatch for umi correction, default 1
  --sat_file TEXT Needs: --umi_on       Output sequencing saturation file, default None
  --scrna,--scRNA,--SCRNA               Set scRNA mode, default false
  --no_filter_matrix Needs: --scrna     Not filter the gene expression matrix, default false

HandleBam version: 1.0.0
```

Required parameters:

* -i filename. Input bam filename
* -o filename. Output bam filename
* -a filename. Input annotation filename
* -s filename. Output summary filename
* -e filename. Output barcode gene expression filename

Optional parameters:

* -q integer. Set mapping quality threshold, default 10
* --save_lq. Save low quality reads(less than paramter of '-q'), default not save
* --save_dup. Save duplicate reads, default not save.
* --anno_mode integer. Select annotation mode, default is 2
* --umi_on. Open umi correction, default off
* --umi_min_num integer. Minimum umi number for correction, default 5
* --umi_mismatch integer. Maximum mismatch for umi correction, default 1
* --sat_file filename. Output sequencing saturation file, depend on --umi_on

### Example

#### Normal Mode

```text
$./install/bin/handleBam \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/exp.summary.txt \
 -e batch25/exp.barcode_gene_exp.txt
```

Input parameters: `-i batch25/batch25.bam -a batch25/batch25.gtf`  
Ouput parameters: `-o batch25/exp.bam -s batch25/exp.summary.txt -e batch25/exp.barcode_gene_exp.txt`

#### Save Mode

```text
$./install/bin/handleBam \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/exp.summary.txt \
 -e batch25/exp.barcode_gene_exp.txt \
 --save_lq \
 --save_dup
```

Save reads with low quality or duplicate, also set bam flags with *BAM_FQCFAIL* or *BAM_FDUP*

#### Set Annotation Mode

```text
$./install/bin/handleBam \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/exp.summary.txt \
 -e batch25/exp.barcode_gene_exp.txt \
 --anno_mode 1
```

Set annotation mode to *Drop-seq V1*, there will be more options to choose from in the future

#### Use Umi Correction

```text
$./install/bin/handleBam \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/exp.summary.txt \
 -e batch25/exp.barcode_gene_exp.txt \
 --umi_on \
 --umi_min_num 10 \
 --umi_mismatch 2
```

Use umi correction for deduplicate

#### Sequencing Saturation

```text
$./install/bin/handleBam \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/exp.summary.txt \
 -e batch25/exp.barcode_gene_exp.txt \
 --umi_on \
 --sat_file batch25/exp.saturation.txt
```

### Result 

The result cantains three files:
* output bam. The bam file contains gene annotation information
* summary. Text file contains some stat metrics of filter, deduplication and annotation,example:

  1. filter and deduplication metrics table

  | TOTAL_READS | PASS_FILTER | UNIQUE_READS | FAIL_FILTER_RATE | DUPLICATION_RATE |
  | ----------- | ----------- | ------------ | ---------------- | ---------------- |
  | 9999404     | 6450176     | 4036230      | 35.4944          | 37.4245          |

  2. annotation metrics table

  Mode 0/1:

  | TOTAL_READS | READS_WRONG_STRAND | READS_RIGHT_STRAND | READ_AMBIGUOUS_GENE_FIXED | AMBIGUOUS_READS_REJECTED |
  | ----------- | ------------------ | ------------------ | ------------------------- | ------------------------ |
  | 4036230     | 194315             | 3841915            | 0                         | 0                        |

  Mode 2:

  | TOTAL_READS | MAP    | EXONIC | INTRONIC | INTERGENIC | TRANSCRIPTOME | ANTISENSE |
  | ----------- | ------ | ------ | -------- | ---------- | ------------- | --------- |
  | 949906      | 949906 | 855200 | 354      | 42033      | 850051        | 9833      |
  | 100.0       | 100.0  | 90.0   | 0.0      | 9.9        | 89.5          | 1.0       |

  3. umi correction metrics table(exists if using umi correction)
  
  umi number stat:

  | BARCODE_GENE_NUM | UMI_CNT_RAW | UMI_CNT_DEDUP | RAW_PCT | DEDUP_PCT |
  | ---------------- | ----------- | ------------- | ------- | --------- |
  | 9055411          | 7972608     | 6861844       | 131.97  | 116.19    |

  umi mismatch positions:

  | POSITION | CNT    | PCT   |
  | -------- | ------ | ----- |
  | 1        | 141568 | 9.67  |
  | 2        | 220451 | 15.06 |
  | 3        | 136248 | 9.31  |
 
  umi mismatch types:

  | TYPE | CNT    | PCT   |
  | ---- | ------ | ----- |
  | A_A  | 0      | 0.00  |
  | A_C  | 108891 | 7.44  |
  | A_G  | 167747 | 11.46 |

  4. sequencing saturation(exists if using parameter '--sat_file')

  | sample | bar_x | bar_y1   | bar_y2 | bin_x | bin_y1   | bin_y2 |
  | ------ | ----- | -------- | ------ | ----- | -------- | ------ |
  | 0.05   | 2     | 0.60219  | 1      | 166   | 0.602246 | 51     |
  | 0.1    | 4     | 0.706105 | 1      | 330   | 0.706156 | 71     |
  | 0.15   | 5     | 0.755016 | 1      | 492   | 0.755062 | 84     |

  sample: sample data size divide by total data size  
  bar_x: using barcode as spot, mean reads per spot  
  bar_y1: sequencing saturation, 1 - (uniq/total)  
  bar_y2: median genes per spot  
  bin_x/bin_y1/bin_y2: using bin number 150 as spot

* gene expression. The text file composed of three columns, which are `Barcode\tGene\tTimes`, example:
 
  | Barcode              | Gene       | Times |
  | -------------------- | ---------- | ----- |
  | TATTGGTACACCTCACCTCC | AL954705.1 | 8     |
  | TCTATCGGGCCTCGCTGTGC | AL954705.1 | 16    |
  | ATGACACCACCGTCTTCCCT | AL954705.1 | 1     |

  Also output matrix market files: *barcodes.tsv.gz genes.tsv.gz matrix.mtx.gz*
