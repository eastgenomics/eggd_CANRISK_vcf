# eggd_CANRISK_vcf

## What does this app do?

This app checks whether variants detected at a set of PRS positions are sufficiently covered and optionally whether they overlap any copy number variants,
which is noted in text output files. Furthermore, the app generates a VCF formatted to be compatible with the [CanRisk tool](https://www.canrisk.org/canrisk_tool/).

## What inputs are required for this app to run?

- sample VCF from Sentieon Haplotyper (called against PRS variants)
- segments VCF from GATK gCNV (optional)
- a read depth threshold can be set (optional, defaults to 20)

## How does this app work?

The app normalises multiallelic variants and modifies the VCF to include ALTs in 0/0 positions.
If required, it also checks the segments file using bedtools intersect, to see if any relevant variants are in CNVs.
Read depth of the called variants is checked and low covered ones are reported in a text output file.


## What does this app output?

- a filtered VCF compatible with CANRISK, excluding variants that do not reach the coverage threshold specified
- a text file listing excluded variants that are covered at less than input depth (default 20x)
- a text file listing PRS positions that are affected by an overlapping CNV

## What limitations does this app have?

This assumes the provided VCF is already limited to the PRS variants list, and contains a genotype for all those variants.

## This app was created by East GLH
