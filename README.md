# Scripts used in the paper Stable Transgenerational Inheritance of DNA methylation

## Set up

* All python scripts are meant for Python 3.4+
* All scripts import sys, math, glob, multiprocessing, subprocess, os, bisect, and random (all/most come up python)
* Many scripts import pandas, numpy, scipy, and/or sklearn
* To avoid package dependencies, use install [anaconda](https://www.continuum.io/downloads)

* Most R packages can be downloaded through CRAN.

## Program versions
Listed below are program versions used for analysis.
* Python: 3.5.2
* Anaconda: 4.1.6
	- Numpy: 1.11.0
	- Scipy: 0.17.1
	- Scikit: 0.17.1
	- Pandas: 0.17.0
* R: 3.2.4
	- ggplot2: 2.2.1
	- reshape2: 1.4.2
	- ply: 1.8.4
	- dplyr: 0.5.0
	- userfriendlyscience: 0.5-2
	- RVAideMemoire: 0.9-62
	- grid 3.2.4
	- gridExtra 2.2.1

## Utility Scripts
### bioFiles.py
* for handling common file types
* needs to be on the python path or in the same directory as other scripts for all/most other scripts to run correctly

### bth_util.py
* extra utility functions
* needs to be on the python path or in the same directory as other scripts for all/most other scripts to run correctly

### combine_allc_pe.py
Combine multiple allC files at basepair level into one allC file

```
Usage: python combine_allc_pe.py [-f] [-p=num_proc] [-o=out_id] [-c=chrm_list | -cf=fasta_index] <allc_path> <sample_name> [sample_name]*
Required:
allc_path    path to allC files
sample_name    name of sample; used to find allC files
               when "-f" flag set, file with sample names listed one per line
Optional:
-f    sample names are in file
-p=num_proc    number of processors to use [default 1]
-o=out_id     output file identifier [default "combined"]
-c=chrm_list    comma-separated list of chrms to use
-cf=fasta_index    fasta index file with chrms to use

```

### filter_allc_coverage_pe.py
Creates new allc files for all input samples that only includes information about positions which have at least minCov reads for each sample

expects all chromosomes in one allC file

```
Usage: python filter_allc_coverage.py [-v=min_cov] <allc_path> <sample1> [sampleN]*
Required:
allc_path    path to allC files
sampleN    name of sample; used to find allC files
Optional:
-v=min_cov    min coverage for positions to include [default 3]
-p=num_proc    number of processors to use [default 1]
```

### unmethylate_allc_pe.py
Creates pseudo-allC file where all positions are unmethylated

Output file has same name as input file with "-unmethylated" appended

```
Usage: python unmethylate_allc_pe.py [-f] [-p=num_proc] [-v=NA] <allc_file> [allc_file]*
Required:
allc_file    allC file to unmethylated
             when "-f" set, file with list of allC files
Optional:
-f    allC files names listed in the file
-p=num_proc    number of processors to use [default 1]
-v=coverage    coverage for each position [default as-is in input]
```

## Transgenerational MA lines

### dmr_gen_counts.py
Used for between-generation computations of methylation

expects all chromosomes in one allC file and when minCov is set, minCov is output from `filter_allc_coverage.py`

```
Usage: python dmr_gen_counts.py [-o=outID] [-m=methType] [-p=numProc] [-v=min_cov]
 <dmrFile> <allcPath> <sample1> <sample2> [sampleN]*
Required:
dmrFile   tab-delimited file (BED format) with DMRs to investigate
allcPath    Path to allC files; all chrms together for each sample
sample    sample names as part of the allC file
Optional:
-o=outID   identifier for output file [default "out"]
-m=methType    methylation type [default C]
-p=numProc    num. of processors to use [default 1]
-v=minCov   min. coverage used as part of allC file name [default None]
```

### dmr_gen_switches.py
Identifies significant methylation changes between generations using Fisher's exact test and minimum change in methylation

Input file is output from `compare_dmrs_gens_pe.py

```
Usage: python dmr_gen_switches.py [-wm] [-n=num_c_thresh] [-m=meth_thresh] [-d=length_thresh] [-f=fdr] [-o=outID] <in_file>
Required:
in_file	tab-delimited file of DMRs and read counts
Optional:
-wm    methylation threshold is for raw methyl difference not percent difference
-n=num_c_thresh    min number of cytosines in region to be considered for analysis [default 10]
-d=lenth_thresh    min length of dmr in bp [default 40]
-m=meth_thresh    min methylation change btwn generations to be considered a switch [default 0.3]
-f=fdr    FDR value for significant switches [default 0.05]
-o=out_id    identifier for output files [default uses input file name]

```

### dmr_file_to_bed.py
Converts the switches output of `dmr_gen_switches.py` to BED file

```
Usage: python dmr_file_to_bed.py [-v=score_thresh] [-p=name_prefix] [-o=outID] <in_file>
Required:
in_file   input file of DMRs
Optional:
-v=score_thresh    min score to include in output [default -1, no threshold]
-p=name_prefix    prefix for naming features [default None]
-o=outID    identifier for output file [default uses input file name]
```

### dmr_counts_pe.py
Computes weighted methylation over regions

```
Usage: python dmr_counts_pe.py [-o=outID] [-m=methTypes] [-p=numProc] [-v=minCov] <dmrFile> <allcPath>  <sample1> <sample2> [sampleN]*
Required:
dmrFile   tab-delimited file (BED format) with DMRs to investigate
allcPath    Path to allC files; all chrms together for each sample
sample    sample names as part of the allC file
Optional:
-o=outID   identifier for output file [default "out"]
-m=methType    methylation context to include [default C]
-p=numProc    number of processors [default 1]
-v=minCov    min. coverage used as part of allC file name [default None]

```

## Epigenotyping

### find_all_mpos_dif_pe.py
get individual positions that differ based on binomial test

Uses allC files specific to each chromosome

```
Usage: python find_all_mpos_dif_pe.py [-v=min_cov] [-c=chrm_list] [-o=out_id] [-p=num_proc] [-m=meth_types] <allc_path> <sample1_name> <sample2_name>
Required:
allc_path    path to allc files
sample_name    names of samples to compare
Optional
-v=min_cov    min coverage to include a position [default 3]
-o=out_id    string for output file name [default "out"]
-c=chrm_list    comma-separated list of chrms [default arabidopsis]
-p=num_proc    num processors to use [default 1]
-m=meth_types    comma-separated list of "CG", "CHG", and/or "CHH" [default all]
```

### filter_pos_gene_gbm.py
Filter a list of positions by gene-body methylation and/or CDS

```
Usage: python filter_pos_gene_gbm.py [-cds] [-v] <pos_file> <gbm_file> <gff_file>
Required:
pos_file    position file, tab-delimited BED format, to be filtered
gbm_file    file with list of gbM genes, one gene per line
            use "none" or "na" to use all genes
gff_file    GFF formatted file with genes
Optional:
-cds    use CDS annotation not gene
-v    include coordinates opposite of what is specified

```