# eggd_CANRISK_vcf

## What does this app do?

This app checks whether variants detected at a set of PRS positions are sufficiently covered and optionally whether they overlap any copy number variants,
which is noted in text output files. Furthermore, the app generates a VCF formatted to be compatible with the [CanRisk tool](https://www.canrisk.org/canrisk_tool/).

## What inputs are required for this app to run?

- sample VCF from Sentieon Haplotyper (called against PRS bed as targets)
- BED of positions relevant for PRS from CanRisk (https://canrisk.atlassian.net/wiki/spaces/FAQS/pages/35979266/What+variants+are+used+in+the+PRS)
- segments VCF from GATK gCNV (optional)

## How does this app work?

The app normalises multiallelic variants and modifies the VCF to include ALTs in 0/0 positions.
If required, it also checks the segments file using bedtools intersect, to see if any relevant variants are in CNVs.


## What does this app output?

- a modified VCF compatible with CANRISK, with "PRS" appended to the file name.
- a text file listing PRS positions that are not covered at 20x
- a text file listing PRS positions that are affected by an overlapping CNV

## What limitations does this app have?

This assumes the provided VCF is already limited to the PRS variants list, and contains a genotype for all those variants.

## This app was created by East GLH
