+++
date = "2017-06-02T00:11:02+01:00"
title = "Epigenotyping"
weight = 40
+++

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
- All of the allC files are named as "allc_samplename_chromosome.tsv".
- All of the allC files are in one directory.

### Determining the correct parameters
- Information coming soon

## Workflow

### Step 1: Find informative cytosines
There should be two parent samples and need to determine which positions have different methylation. Assuming you have allC files, this is determine by comparing the final column of the allC file which says if the cytosine is methylated and positions which is methylated in one sample but not in the other are included.

Note that low-covered positions are considered "unmethylated", so the coverage parameter is important to make sure positions are differentially methylated because there is differential methylation not due to lack of coverage.

In the *sample_data* folder, run [find_mpos_diff_pe.py](/appendix/#find-mpos-diff-pe-py)

```bash
python ../scripts/find_mpos_diff_pe.py -v=3 -c=Chr1 -o=sample-data allc_files \
mother-sample father-sample
```

This will create a file named *sample-data_mpos_diff.tsv*.

**Note**: If methylation data does not include the final column which determines if site is methylated, use a different criteria to determine which sites are differentially methylated between the samples. Put the results in a tab-delimited file with two columns: chromosome and position. This file can use used in the next step instead of the output from `find_mpos_diff_pe.py`.

### Step 1.5: Filter by feature type [optional]

For samples such as epiRILs, it's important to filter cytosines by feature type, such as genes and/or gene-body methylated genes. Lists of gene-body methyalted genes can be found from [Bewick et al. (2016)](http://www.pnas.org/content/113/32/9111.long). Positions could be filtered by other characteristics as well.

The script requires the position list to filter and a GFF file with gene annotation. It can also accept a gene subset list. If you don't want to subset the genes, specificy `na` instead of subset file. By default, positions within the genes are included. Use `-cds` to only have positions in CDS regions. 

The `-v` option selects opposite of genes/subset. If no gene subset is given `-v` gets all positions outside of genes/CDS. If gene subset is given, `-v` gets all positions of genes not in the subset list.

Data in not included but could be run. In the *sample_data* folder, run [filter_pos_gene_subset.py](/appendix/#filter-pos-gene-subset-py)

```bash
python ../scripts/filter_pos_gene_subset.py sample-data_mpos_diff.tsv <subset_file> <gff_file>
```

Output (*sample-data_mpos_diff_subset_gene.tsv*) would have all positions in of entire gene body of genes listed in subset file.

```bash
python ../scripts/filter_pos_gene_subset.py -v sample-data_mpos_diff.tsv na <gff_file>
```
Output (*sample-data_mpos_diff_nongenes_gene.tsv*) would have all non-genic positions.

### Step 2: Methylation level of samples
Once we have the informative positions, we need the methylation level of all samples we want to analyze. Again, the coverage parameter is important to get an accurate methylation level at each position.

In the *sample_data* folder, run [weighted_meth_by_pos_pe.py](/appendix/#weighted-meth-by-pos-pe-py)

```bash
python ../scripts/weighted_meth_by_pos_pe.py -v=3 sample-data_mpos_diff.tsv \
allc_files mother-sample father-sample F2-1 F2-2 F2-3 F2-4 F2-5
```

This will create a file named *sample-data_wm_pos_Chr1.tsv*.

**Note**: allC files used by this script do not need to contain the final "isMethylated" column. They must follow the same naming convention and have the following columns: chromosome, position, strand, 3-bp sequence context, methylated read count, total read count.

### Step 3: Epigenotyping
The main input for epigenotyping is the output of `weighted_meth_by_pos_pe.py`.

This is the main script. It accepts several parameters. Defaults were chosen for the paper but will need to be updated depending on species, filial generation, coverage, ect. 

In the *sample_data* folder, run [epigenotyping_pe_v7.3.py](/appendix/#epigenotyping-pe-v7-3-py)

```bash
python ../scripts/epigenotyping_pe_v7.3.py -g=2 -b=50kb -c=3 -t=3300kb,6300kb \
-m=mother-sample -f=father-sample sample-data_wm_pos_Chr1.tsv
```

This will create a file name *sample-data_epigenotype-v7.3_50kbp_g-2_Chr1_cb-3_both.tsv*

### Step 4: Graphing
Information coming soon.

### Step 5: Finding crossovers
Once you have an epigenotype map, you can get approximate crossover locations and number of crossovers.

In the *sample_data* folder, run [find_crossovers.py](/appendix/#find-crossovers-py)

```bash
python ../scripts/find_crossovers.py sample-data_epigenotype-v7.3_50kbp_g-2_Chr1_cb-3_both.tsv
```

This will create two files: *sample-data_epigenotype-v7.3_50kbp_g-2_Chr1_cb-3_both_count.tsv* (number of breakpoints per sample) and *sample-data_epigenotype-v7.3_50kbp_g-2_Chr1_cb-3_both_pos.tsv* (location of breakpoints)

## SNP Verification


If SNP information is available, it is possible to use WGBS to to verify the epigenotype map. 

### Step 1: Samtools mPileup

This is require a BAM file of the mapped reads and a BED formatted file of SNP locations then running [Samtools](http://www.htslib.org/doc/samtools-1.2.html) mpileup getting the reads at the SNP locations. Remember, BED files are 0-based and half-open (start position is included, end position not). If the SNP occurs at 1-based position *x*, the BEC file coordinates would be *x-1, x*.

```bash
samtools mpileup -l bed_file bam_file > out.txt
```

BAM files are not included in the sample data due to file size. Pileup files are available in *sample_data/pileup_files*.

### Step 2: Decode pileup files

Samtools pileup files use non-obvious notation for the read-level nucleotides, so this script decodes the output into a user-friendly format. 

If interested, `.` stands for a match to the reference on the forward strand, `,` for a match on the reverse strand, `ACGTN` for mismatch on forward strand, and `acgtn` for mismatch on reverse strand.

In *sample_data*, run [decode_pileup_pe.py](/appendix/#decode-pileup-pe-py),

```bash
python ../scripts/decode_pileup_pe.py -o=snp-info pileup_files/F2-1_pileup.txt \
pileup_files/F2-2_pileup.txt pileup_files/F2-3_pileup.txt pileup_files/F2-4_pileup.txt \
pileup_files/F2-5_pileup.txt pileup_files/mother-sample_pileup.txt pileup_files/father-sample_pileup.txt
```

This will output *snp-info_decoded_pileup.tsv*.

The output file includes sample name, chromosome, position, and 8 columns "A,C,G,T,a,c,t,g". "A,C,G,T" are forward strand reads and "a,c,g,t" are reverse strand reads. Number indicates number of reads supporting that nucleotide at the position.

### Step 3: Determine genotype

With read information from all SNP positions, we can how genotype the samples. It determines which nucleotides are different between the mother and father sample and guesses genotype. Positions where mother and father are inditinguishable are reported as such.

For example, if the mother sample has reads supporting "A" and "c" and the father sample has reads supporting "A" and "g", "c" and "g" will be used to determine genotype. Samples with only "c" reads will be called "mother", samples with only "g" reads will be called "father", and samples with "c" and "g" reads will be called "heterozygous".

In the *sample_data* folder, run [pileup_genotype_pe.py](/appendix/#pileup-genotype-pe-py),

```bash
python ../scripts/pileup_genotype_pe.py -m=mother-sample -f=father-sample \
snp-info_decoded_pileup.tsv
```

This creates output file *snp-info_decode_pileup_genotyped.tsv*.

**Note**: In the output file, do not use positions with have "NA" listed as the genotype for any of the samples. 
When one parent is genotyped "NA", there was a distinguishable difference between the parents but one parent was incorrectly genotype (this can happen when one parent is heterozygous). 
Positions where one of the other samples is genotype "NA" indicates insufficient coverage at that position to determine genotype.

For the sample data, the first two SNPs are indistinguishable in the parents. The last SNP has multiple "NA" genotypes. The third SNP is the only useful position.