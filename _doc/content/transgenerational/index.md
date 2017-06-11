+++
date = "2017-06-02T00:11:02+01:00"
title = "Transgenerational Inheritance"
weight = 30
+++

**Sample data is not provided for this topic.**
For illustrative purposes, assume that three generations of line "A" are available and they are named "A-1", "A-2", "A-3". This species only has two chromosomes, "Chr1" and "Chr2".

## Beginnings

### Set up requirements
- All samples have allC files.
    - AllC file are tab-delimited.
    - They have the following columns
        - Chromosome
        - Position (1-indexed)
        - Strand (+ or -)
        - Context (3 bp sequence)
        - Number methylated reads
        - Number total reads
        - Is the position methylated? (1 for yes)
- All of the allC files are named as "allc_samplename_chromosome.tsv" or "allc_samplename.tsv".
- Many of the following scripts expect a single allC file per sample. In this case, the file is named "allc_samplename.tsv". Create this file using `cat allc_samplename_*.tsv > allc_samplename.tsv`.
- All of the allC files are in one directory.

## Identification of Methylated Regions

### Step 1: Combine samples allC files
We combine the information from all samples in a line to form a "pan-methylome". To do this, we sum up the number of methylated and number of total reads at each position across all samples and write the results to a new allC file.

For example run, [combine_allc_pe.py](appendix/#combine-allc-pe-py)

```
    python combine_allc_pe.py -o=lineA-pan-methylome -c=Chr1,Chr2 allc_path A-1 A-2 A-3
```

This creates two files (one for each chromosome) *allc_lineA-pan-methylome_Chr1.tsv* and *allc_lineA-pan-methylome_Chr2.tsv*.

### Step 2: Unmethylate the pan-methylome
We are trying to identify all methylated regions, so we want to create a null, unmethylated genome to compare it against.

```
    python unmethylate_allc_pe.py allc_lineA-pan-methylome_Chr1.tsv allc_lineA-pan-methylome_Chr2.tsv
```

This creates two files (one for each input): *allc_lineA-pan-methylome-unmethylated_Chr1.tsv* and *allc_lineA-pan-methylome-unmethylated_Chr2.tsv*

### Step 3: DMR analysis
Using your prefer DMR analysis tool, you can compare the pan-methylome and unmethylated pan-methylome to identify methylated regions of the genome. Parameters and filters should be the same for this analysis as used for between-samples DMR analysis.

## Epiallele Identification

### Step 1: Putative DMR list

Using the preferred DMR analysis tool, identify DMRs between all samples in a line. Depending on the program, this DMR list should be filtered for length, methylation change, ect. 
After filtering, we have a set of putative epialleles. For illustrative purposes, these regions are saved in a file named *dmr_list.tsv*.

### Step 2: Coverage filters

Many of the scripts in this analysis accept an optional minimum coverage parameter which only keeps positions with sufficient coverage in all samples. We apply this minimum coverage to minimize methylation variation of a region/sample due to missing information.

To save computational time later, we will create a new allC file for each sample which only includes positions with sufficient coverage. This eliminates needing to apply this filter repetitively.

For example, run [filter_allc_coverage_pe.py](appendix/#filter-allc-coverage-pe-py) to with minimum coverage 2. This script expects all chromosomes in one allC file.

```
    python filter_allc_coverage_pe.py -v=2 allc_path A-1 A-2 A-3
```

This creates files *allc_A-1_cov2.tsv*, *allc_A-2_cov2.tsv*, and *allc_A-3_cov2.tsv*.

### Step 3: Between generation comparison

With the DMR list, we want to get the number of methylated and unmethylated reads for each region and compare adjacent generations. First, we will get the read counts for organized by each comparison. When specifying samples, **order matters**.

To allow for fair comparsion between samples, we want to use positions which have sufficient coverage in all samples. So we use the coverage allC files created in Step 2. In this script, the coverage parameter is used to search for files named as `allc_samplename_cov#.tsv`, which includes all chromosomes.

Run [dmr_gen_counts.py](appendix/#dmr-gen-counts-py) for minimum coverage 2.

```
    python dmr_gen_counts.py -v=2 dmr_list.tsv allc_path A-1 A-2 A-3
```

This creates file **.

### Step 4. Identify signficant DMRs
