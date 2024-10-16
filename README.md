# eggd_CANRISK_vcf

## What does this app do?

First, it checks for PRS variants which were not able to be called (i.e. that have a GT = ./.) and records these variants in a .txt file. If no variants are found no .txt is outputted. If requested via the `convert_gt_no_call` input (default is `true`), their genotypes are converted to 0/0.

Next, read depths of PRS variants are inspected and variants below the threshold are reported in a text output file. If no PRS variants with sub-threshold read depth exist, then no coverage check text file is outputted. Optionally, these variants with sub-threshold read depth may have their genotype converted to 0/0, as specified by the `convert_gt_low_dp` input (default is `true`).

After this, it checks the CNV segments file (if provided) to see if any PRS variants are in CNVs. Details of these variants are outputted in a text file, if no such variants are present then no text file is outputted. PRS variants which occur in CNVs can optionally have their genotype converted to 0/0 as specified by the `convert_gt_cnv` input (default is `true`).

Lastly, the app outputs a VCF formatted to be compatible with the [CanRisk tool](https://www.canrisk.org/canrisk_tool/).

## What inputs are required for this app to run?

- sample VCF from Sentieon Haplotyper (called against PRS variants)
- segments VCF from GATK gCNV (optional)
- a read depth threshold can be set (optional, defaults to 20)
- boolean whether to convert the genotype of PRS variants which intersect with CNVs to 0/0 (optional, default true)
- boolean whether to convert the genotype of PRS variants with low coverage (as specified by the depth threshold input) to 0/0 (optional, default true)
- boolean whether to convert the genotype of PRS variants that were not called (i.e. GT = ./.) to 0/0 (optional, default true)


## How does this app work?
The app first identifies variants which were unable to be called via `bcftools filter` using the condition `FORMAT/GT=="./."`. If such variants are found, a .txt file containing a description of the variants (CHROM, POS, REF, ALT, DP, GT) is outputted using `bcftools isec` and `bcftools query`. If requested, the genotypes of these variants in the sample VCF are converted to 0/0 using a custom bash function `convert_gt` (see eggd_CANRISK_vcf.sh for more details).

Next, the script finds variants in the sample VCF with sub-threshold depth values using `bcftools filter` and the specified depth threshold given via the `depth` input. If requested, the app will then convert the GT values for these variants to 0/0 using the `convert_gt` function and a .txt file describing any such variants is outputted (as described above).

If a segments VCF is provided via the `segments_vcf` input, it checks for any CNVs that have been called via `bcftools filter` and then for any overlaps of the CNVs with PRS variants in the sample VCF using `bcftools isec`. If requested, it converts the GTs of CNV overlapping variants to 0/0 via the `convert_gt` function and a .txt file is outputted describing these variants (as described above).

## What does this app output?

- a filtered VCF compatible with CANRISK
- a text file listing variants that are covered at less than input depth (default 20x), if no such variants exist then no text file is outputted.
- a text file listing PRS variants that are affected by an overlapping CNV, if no such variants exist then no text file is outputted.
- a text file listing PRS variants that were not able to be called (i.e. their GT = ./.), if no such variants exist then no text file is outputted.

## What limitations does this app have?

This assumes the provided VCF is already limited to the PRS variants list.

## This app was created by East GLH
