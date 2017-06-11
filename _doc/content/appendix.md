+++
date = "2017-06-02T00:11:02+01:00"
title = "Appendix"
weight = 70
+++

This page includes list of inputs and parameter options for all scripts. This should not be used as workflow or tutorial, but is a resource when necessary for specific scripts.

## combine_allc_pe.py
- Require imports: sys, multiprocessing, subprocess, os, gzip
- Expects allC files to be for a single chromosome
- Accepts gzipped allC files

```
Usage:  python combine_allc_pe.py [-q] [-h] [-f] [-p=num_proc] [-o=out_id] 
        [-c=chrm_list | -cf=fasta_index] <allc_path> <sample_name> [sample_name]*

Required:
allc_path       path to allC files
sample_name	    name of sample; used to find allC files
                when "-f" flag set, file with sample names listed one per line

Optional:
-q              quiet; do not print progress
-h		        print help message and exit
-f		        sample names are in file
-p=num_proc	    number of processors to use [default 1]
-o=out_id	    output file identifier [default "combined"]
-c=chrm_list	comma-separated list of chrms to use [default for arabidopsis]
-cf=fasta_index	fasta index file with chrms to use
```

## unmethylate_allc_pe.py
- Required imports: sys, multiprocessing, subprocess, os, gzip
- Accepts gzipped allC files
- **Note**: specify allC file not allC path and sample names
- Output files are named as input file name with `unmethylated` added

```
Usage:	python unmethylate_allc_pe.py [-q] [-h] [-f] [-p=num_proc] [-v=coverage]
	    <allc_file> [allc_file]*

Required:
allc_file     allC file to unmethylated
              when "-f" set, file with list of allC files

Optional:
-q            quiet; do not print progress
-h		      print help message and exit
-f 		      allC files names listed in the file
-p=num_proc   number of processors to use [default 1]
-v=coverage   coverage for each position [default as-is in input]
```

## filter_allc_coverage_pe.py
- Required imports: sys, math, glob, multiprocessing, subprocess, os
- Other required files: bioFiles
- Expects allC files to have all chromosomes for one sample

```
Usage:  python filter_allc_coverage_pe.py [-q] [-h] [-v=min_cov] <allc_path>
	    <sample_name> [sample_name]*

Required:
allc_path	 path to allC files
sampleN		 name of sample; used to find allC files

Optional:
-q		     quiet; do not print progress
-h		     print help and exit
-v=min_cov   min coverage for positions to include [default 3]
-p=num_proc  number of processors to use [default 1]
```



## find_mpos_diff_pe.py
 - Required imports: sys, multiprocessing, subprocess, os, gzip
 - Expects allC files to be split by chromosomes
 - Accepts gzipped allC files

```
Usage:	python find_mpos_diff_pe.py [-h] [-q] [-v=min_cov] [-c=chrm_list]
    [-o=out_id] [-p=num_proc] [-m=meth_types] <allc_path> <sample1_name>
	<sample2_name>

Required:
allc_path	path to allc files
sample_name	names of samples to compare

Optional:
-h      print this help screen and exit
-q      quiet; does not print progress
-v=min_cov      min coverage to include a position [default 3]
-o=out_id       string for output file name [default "out"]
-c=chrm_list    comma-separated list of chrms [default arabidopsis]
-p=num_proc     num processors to use [default 1]
-m=meth_types   comma-separated list of "CG", "CHG", and/or "CHH"
                [default all]

```

## filter_pos_gene_subset.py
- Required imports: sys, os, bisect
- Input should be the output of *find_mpos_diff_pe.py* to filter positions

```
Usage:	python filter_pos_gene_subset.py [-h] [-q] [-cds] [-v] <pos_file>
	    <subset_file> <gff_file>

Required:
pos_file	 1-based position file, tab-delimited columns: chrm, start, end
subset_file	 file with list of genes to subset, one gene per line
		     use "none" or "na" to use all genes
gff_file	 GFF formatted file with genes

Optional:
-h		     print help and exit
-q		     quiet; do not print progress
-cds		 use CDS annotation not gene
-v		     include coordinates opposite of what is specified
```

## weighted_meth_by_pos_pe.py
- Required imports: sys, multiprocessing, subprocess, os, bisect, gzip
- Accepts gzipped allC files

```
Usage:	python find_mpos_diff_pe.py [-h] [-q] [-v=min_cov] [-c=chrm_list]
	    [-o=out_id] [-p=num_proc] [-m=meth_types] <allc_path> <sample1_name>
	    <sample2_name>

Required:
allc_path	path to allc files
sample_name	names of samples to compare

Optional:
-h		print this help screen and exit
-q		quiet; does not print progress
-v=min_cov	min coverage to include a position [default 3]
-o=out_id	string for output file name [default "out"]
-c=chrm_list	comma-separated list of chrms [default arabidopsis]
-p=num_proc	num processors to use [default {:d}] 1
-m=meth_types	comma-separated list of "CG", "CHG", and/or "CHH"
		[default all]
```

## epigenotyping_pe_v7.3.py
- Required python imports: sys, math, multiprocessing, subprocess, os, numpy, pandas, sklearn, functools
- Other required files: bth_util, decodingpath, transitions
    - decodingpath.py includes code for the smoothing algorithms (forward-backward and Viterbi)
    - transitions includes code to compute the transition matrix needed by decoding algorithms
- This script is run *per chromosome*.

```
Usage:	python epigenotyping_pe_v7.3.py [-q] [-n-mpv] [-t-out] [-g=generation]
	[-c=bin_thresh] [-d=decoding_type] [-p=num_proc] [-o=out_id] [-m=mother_
	samples][-f=father_samples] [-b=bin_size] [-t=centromere] <input_file>

Requried:
input_file	tab-delimited file of of weighted methylation by position for samples

Optional:
-q		quiet; do not print progress
-n-mpv		do not check for systematic mid-parent bias
-t-out		write transition matrix to file
-g=generation	generation of self-crossing; used to determine classification
		probabilities; use 0 for uniform weight [default 2]
-d=decode_type	decoding type to use (capitlization ignored) [default A]
		Viterbi="v" or "viterbi"
		Forward-Backward="forwardbackward", "f" or "fb"
		Both="all" or "a"
		Off="false", "none", or "n"
-o=out_id	identifier for output file [default "out" or variation of
		input file name]
-p=num_proc	number of processors [default 1
-c=bin_thresh	minimum number of features per bin to be classified
		groups bins to reach this number [default 3
-m=mother_samples	comma-separated sample name(s) of mother
		[default mother]
-f=father_samples	comma-separated sample name(s) of father
		[default father]
-b=bin_size	size of bins in bp [default 100kbp]
-t=centromere	centromere coordinates as "start,end"; can include multipe
		centromeres as "start1,end1,start2,end2..." [default None]
```

## find_crossovers.py
- Required imports: sys, os, pandas
- This script is run per chromosome.

```
Usage:	python find_crossovers.py [-c=predict_column] [-o=out_id] <input_file>

Required:
input_file	tab-delimited file with samples epigenotype per bin

Optional:
-h		print this help menu and exit
-q		quiet;do not print progress
-o=out_id	output identifier [default variation of
		input file name
-c=predict_column 	label of column to use as final epigenotype
		[default "vit.prediction"]
```

## decode_pileup_pe.py
- Required imports: sys, multiprocessing, subprocess, os
- Assumes all pileup files have information for the same positions.

``` 
Usage: python decode_pileup_pe.py [-o=out_id] [-p=num_proc] <pileup_file> [pileup_file]*

Required:
pileup_file	pileup file for a sample; output from samtools pileup

Optional:
-q		quiet; do not print progress
-o=out_id	identifier for output file [default "out"]
-p=num_proc	number of processors [default 1]

```

## pileup_genotype_pe.py
- Required imports: sys, multiprocessing, subprocess, os
- Uses the output of *decode_pileup_pe.py*

```
Usage:	python pileup_genotype_pe.py [-h] [-q] [-o=out_id] [-p=num_proc]
	[-m=mother_label] [-f=father_label] [-v-min_cov] <decoded_pileup_file> 

Required:
decode_pileup_file	tab-delimited input file; output of decode_pileup_pe.py

Optional:
-h		print help and exit
-q		quiet; do not print progress
-o=out_id	identifier for output file [default variation of input file
		name]
-v=min_cov	min number of reads needed to support genotype [default 1]
-p=num_pro	number of processors [default 1]
-m=mother_label	sample name of mother [default mother]
-f=father_label	sample name of father [default father]
```