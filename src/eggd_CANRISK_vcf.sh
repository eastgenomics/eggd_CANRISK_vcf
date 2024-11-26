#!/bin/bash

set -exo pipefail

convert_gt () {
    # Convert the genotype for variants in the given regions file to 0/0.
    original_vcf=$1
    regions_file=$2

    # create intersect vcf using the provided region and vcf files
    bcftools isec "$original_vcf" -T "$regions_file" -w1 | \
        bcftools view -H -o intersect.tsv

    # make vcf containing non-selected variants, the ^ causes everything which
    # does not intersect the regions file to be used
    bcftools isec "$original_vcf" -T ^"$regions_file" -w1 -o intersect_removed.vcf

    # convert selected variant's genotypes to 0/0
    paste <(cut -f1-9 intersect.tsv) \
        <(cut -f10 intersect.tsv | sed 's/^[^:]*/0\/0/') \
        <(cut -f11- intersect.tsv) > converted_intersect.tsv

    # add the converted variant entries to the non-selected variants vcf and sort
    # intermediate "to_be_checked.vcf.gz" file created so as to not overwrite the original
    # vcf file until checks have been passed
    cat intersect_removed.vcf converted_intersect.tsv | bcftools sort -o to_be_checked.vcf.gz
    tabix -f to_be_checked.vcf.gz

    # check that all other columns have remained unmodified
    if ! cmp -s <(bcftools view "$original_vcf" -H | cut -f1-9,11-) \
        <(bcftools view to_be_checked.vcf.gz -H | cut -f1-9,11-); then
            echo "ERROR: conversion of GT values affected integrity of VCF"
            exit 1
    fi

    # check that all GTs have been converted
    gt=$(bcftools isec to_be_checked.vcf.gz -T "$regions_file" -w1 | \
        bcftools query -f '[%GT\n]' | sort | uniq)
    if [[ "$gt" != "0/0" && "$gt" != "0" ]]; then
        echo "ERROR: failed to convert GT values"
        exit 1
    fi

    # now that output file has been checked overwrite/update original file and
    # re-index original vcf
    mv to_be_checked.vcf.gz "$original_vcf"
    tabix -f "$original_vcf"

    # clean up intermediate files
    rm intersect.tsv intersect_removed.vcf converted_intersect.tsv to_be_checked.vcf.gz.tbi
}

main() {

    ### INPUTS
    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    ### PROCESS
    # check for input conflicts
    if [[ -z "$segments_vcf_path" ]] && $convert_gt_cnv; then
        echo "ERROR: no segments file was provided via segments_vcf but convert_gt_cnv was set to true"
        exit 1
    fi

    # get full sample name
    sample_name=$(basename "$sample_vcf_path" | awk -F '_' '{ print $1 }')
    tabix -f "$sample_vcf_path"

    mark-section "Checking for uncalled PRS variants"
    bcftools filter -i 'FORMAT/GT=="./."' "$sample_vcf_path" | \
        bcftools query -f '%CHROM\t%POS\n' > uncalled_coords.tsv

    if [ -s uncalled_coords.tsv ]; then
        mark-section "Writing uncalled variant check file"
        echo "The following PRS variants were not called in this sample:" > "$sample_name"_uncalled_variants_check.txt
        echo -e "\n#CHROM\tPOS\tREF\tALT\tGT" >> "$sample_name"_uncalled_variants_check.txt
        bcftools isec "$sample_vcf_path" -T uncalled_coords.tsv -w1 | \
                bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' >> "$sample_name"_uncalled_variants_check.txt

    else
        echo "No uncalled variants found in the sample VCF"
    fi

    mark-section "Checking for low coverage"
    # identify variants with low coverage (below depth threshold input)
    filter="FORMAT/DP<$depth"
    bcftools filter -i "$filter" "$sample_vcf_path" | \
        bcftools query -f '%CHROM\t%POS\n' > low_dp_coords.tsv

    if [ -s low_dp_coords.tsv ]; then

        mark-section "Writing coverage check file"
        # write list of affected PRS variants to output file
        bcftools filter -i "$filter" "$sample_vcf_path" > low_cov_PRS.vcf
        echo "The following PRS variants are not covered to $depth x read depth:" > "$sample_name"_coverage_check.txt
        echo -e "\n#CHROM\tPOS\tREF\tALT\tDP\tGT" >> "$sample_name"_coverage_check.txt
        # write relevant info about affected variants from filtered VCF
        bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%DP\t%GT\n]' low_cov_PRS.vcf >> "$sample_name"_coverage_check.txt
        echo -e "\nEnd of file" >> "$sample_name"_coverage_check.txt

    else
        echo "No low coverage variants present in sample VCF"
    fi

    if [[ -n "$segments_vcf_path" ]]; then
        mark-section "Checking CNV segments filename matches sample from input VCF"
        segment_sample_name=$(basename "$segments_vcf_path" | awk -F '_' '{ print $1 }')
        if ! [ "$sample_name" == "$segment_sample_name" ]; then
            echo "ERROR: sample names do not match between sample VCF and segment VCF"
            exit 1
        fi

        mark-section "Checking for overlapping CNVs"
        # if segments file input is provided, check if any PRS variants are in a CNV
        # by converting CNV coordinates to a bed file
        # identify CNVs from the segments VCF (bcftools filter)
        # parse their coordinates to a bed file (bcftools query)
        # check for overlaps between that and the sample (PRS) vcf (bcftools isec)
        bcftools filter -e 'FORMAT/GT=="0/0"' "$segments_vcf_path" | \
            bcftools query -f '%CHROM\t%POS\t%INFO/END\n' > CNV_coords.bed
        bcftools isec "$sample_vcf_path" -T CNV_coords.bed -w1 | grep -v ^# > CNV_overlaps.tsv

        if [ -s CNV_overlaps.tsv ]; then
            mark-section "Writing CNV check file"

            # write list of affected PRS variants to output file
            echo "The following PRS variants are located within a CNV call in this sample:" > "$sample_name"_cnv_check.txt
            echo -e "\n#CHROM\tPOS\tREF\tALT\tGT" >> "$sample_name"_cnv_check.txt

            # write relevant info about affected variants
            bcftools isec "$sample_vcf_path" -T CNV_coords.bed -w1 | \
                bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' >> "$sample_name"_cnv_check.txt

            echo -e "\nEnd of file" >> "$sample_name"_cnv_check.txt

        else
            echo "No CNVs called in PRS regions"
        fi
    else
        echo "No segment VCF provided, therefore no CNV checking was performed for this sample."
    fi

    # Modify VCF
    mark-section "Modifying GTs where neccessary"
    if [ -s uncalled_coords.tsv ]; then
        if $convert_gt_no_call; then
            mark-section "Converting genotype of uncalled PRS variants to 0/0"
            convert_gt "$sample_vcf_path" uncalled_coords.tsv
        else
            echo "Genotypes of uncalled PRS variants were not converted to 0/0"
        fi
    fi
    if [ -s low_dp_coords.tsv ]; then
        if $convert_gt_low_dp; then
            mark-section "Converting low coverage PRS variant genotypes to 0/0"
            convert_gt "$sample_vcf_path" low_dp_coords.tsv
        else
            echo "Genotypes of low coverage PRS variants were not converted to 0/0"
        fi
    fi
    if [ -s CNV_overlaps.tsv ]; then
        if $convert_gt_cnv; then
            mark-section "Converting genotype of CNV intersected PRS variants to 0/0"
            convert_gt "$sample_vcf_path" "CNV_coords.bed"
        else
            echo "Genotypes of CNV intersecting PRS variants were not converted to 0/0"
        fi
    fi

    gunzip -c "$sample_vcf_path" > "$sample_name"_canrisk_PRS.vcf

    # Modify variant representation for specific variants to represent them
    # non-parsimoniously as is required by Canrisk
    sed -i 's/4\t84370124\trs10718573\tTA\tT/4\t84370124\trs10718573\tTAA\tTA/' "$sample_name"_canrisk_PRS.vcf
    sed -i 's/22\t38583315\trs138179519\tA\tAAAAG/22\t38583315\trs138179519\tAAAAG\tAAAAGAAAG/' "$sample_name"_canrisk_PRS.vcf

    ### OUTPUTS
    mark-section "Uploading output files"
    # upload files
    canrisk_VCF=$(dx upload "$sample_name"_canrisk_PRS.vcf --brief)
    dx-jobutil-add-output canrisk_PRS "$canrisk_VCF" --class=file

    if [ -s "$sample_name"_uncalled_variants_check.txt ]; then
        uncalled_check=$(dx upload "$sample_name"_uncalled_variants_check.txt --brief)
        dx-jobutil-add-output uncalled_check "$uncalled_check" --class=file
    fi

    if [ -s "$sample_name"_coverage_check.txt ]; then
        coverage_check=$(dx upload "$sample_name"_coverage_check.txt --brief)
        dx-jobutil-add-output coverage_check "$coverage_check" --class=file
    fi

    if [ -s "$sample_name"_cnv_check.txt ]; then
        cnv_check=$(dx upload "$sample_name"_cnv_check.txt --brief)
        dx-jobutil-add-output cnv_check "$cnv_check" --class=file
    fi

	dx-upload-all-outputs --parallel
	mark-success
}
