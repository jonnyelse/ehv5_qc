ehv5_qc

ehv5_qc is a DNAnexus applet for the UK Biobank Research Analysis Platform (RAP) that runs ExpansionHunter v5 on CRAM files and produces annotated, indexed STR genotype VCFs.

This repository contains the implementation used for STR genotyping and quality-control preparation in my PhD thesis. 

Overview

The applet takes one or more aligned CRAM files, a reference genome, and an ExpansionHunter variant catalogue. For each input sample, it:

downloads the CRAM and CRAI index from DNAnexus;
runs ExpansionHunter v5;
sorts and indexes the realigned BAM output;
compresses and indexes the sample VCF;
reheaders the VCF using the requested sample/output prefix;
identifies multi-STR catalogue entries;
adds additional annotation tags used in downstream QC;
uploads the annotated VCF and index to the requested DNAnexus output folder.

The main output is an annotated, bgzipped, tabix-indexed VCF for each analysed sample.

Intended use

This applet was designed for chromosome-scale STR genotyping workflows on the UK Biobank RAP, where ExpansionHunter v5 is run across many samples using a shared variant catalogue.

It is intended as a genotyping and annotation step before downstream sample-level and locus-level STR quality control. 

Inputs
Name	Type	Required	Description
reads	array	Yes	Input CRAM files. Each CRAM is processed as one sample.
reference	file	Yes	Reference genome FASTA used by ExpansionHunter.
variant_catalog	file	Yes	ExpansionHunter variant catalogue.
output_prefix	array	Yes	Output/sample prefix for each input CRAM. These values are also used when reheadering the VCF.
output_folder	string	Yes	DNAnexus folder where output VCFs and indexes are uploaded.
sex	array	No	Sample sex value passed to ExpansionHunter.
threads	int	No	Number of threads passed to ExpansionHunter.
analysis_mode	string	No	ExpansionHunter analysis mode.



Outputs
Name	Type	Description
vcf	array	Annotated, bgzipped ExpansionHunter VCF files.
index_vcf	array	Tabix indexes for the annotated VCF files.


Dependencies

ExpansionHunter v5
bcftools
tabix
samtools
sambamba
jq
R
R packages including dplyr, optigrab, purrr, jsonlite, stringr, tidyr, and data.table
Example DNAnexus usage
dx run ehv5_qc \
  -i reads=sample.cram \
  -i reference=reference.fa \
  -i variant_catalog=variant_catalog.json \
  -i output_prefix=sample_id \
  -i output_folder=/path/to/output \
  -i sex=male \
  -i threads=16 \
  -i analysis_mode=seeking

For multiple samples, provide matching arrays for reads and output_prefix.

Notes on annotation

The applet creates a list of variants from the ExpansionHunter VCF and identifies catalogue entries containing more than one STR motif or repeat structure. These are passed to add_extra_tags.r, which adds additional annotations used by downstream STR QC scripts.

This step was included to support later filtering and interpretation of ExpansionHunter output in large-scale STR analyses.

Scope and limitations

This applet performs STR genotyping and VCF annotation only. It assumes that:

input CRAM files and CRAI indexes are available on DNAnexus;
the reference genome matches the alignment and variant catalogue;
downstream sample-level and locus-level QC will be performed separately;
downstream association testing will be performed by separate applets or scripts.



Reproducibility

This repository contains the source code used to generate ExpansionHunter v5 STR genotype outputs for the analyses described in my PhD thesis. The implementation has been retained to document the computational workflow used at the time of analysis.

Author

Jonny Else
