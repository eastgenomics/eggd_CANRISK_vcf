# eggd_CANRISK_vcf

## What does this app do?

This app formats VCFs to calculate PRS scores using the online CanRisk tool (https://www.canrisk.org/canrisk_tool/).
It also informs the user if any PRS variants are located in a CNV in that sample (optional).

## What inputs are required for this app to run?

- VCF from Sentieon (called against PRS bed as targets)
- PRS file from CanRisk (https://canrisk.atlassian.net/wiki/spaces/FAQS/pages/35979266/What+variants+are+used+in+the+PRS)
- Segments VCF from GATK gCNV (optional)

## How does this app work?

The app modifies the VCF to include ALTs in REF/REF positions & normalises multiallelic variants.
If required, it also checks the intervals file using bedtools intersect, to see if any relevant variants are in CNVs.


## What does this app output?

A single VCF with the PRS name appended to the VCF name.

## What limitations does this app have?

This assumes the provided VCF is already limited to the PRS variants list, and contains a genotype for all those variants.

## This app was created by East Genomics GLH
