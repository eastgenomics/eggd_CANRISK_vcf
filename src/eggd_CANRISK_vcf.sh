#!/bin/bash

set -exo pipefail

main() {
    ### INPUTS

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    ### PROCESS

    # get full sample name
    sample_name=$(basename $sample_vcf_path | awk -F '_' '{ print $1 }')

    ## 1. check CNV overlap:
    # if segments file input is provided, check if any PRS positions are in a CNV
    # by converting CNV coordinates to a bed file
    if [ $segments_vcf_path ]; then
        mark-section "Checking for overlapping CNVs"
        # check CNV VCF filename matches sample from input VCF
        segment_sample_name=$(basename $segments_vcf_path | awk -F '-' '{ print $1"-"$2 }')
        # TODO future proofing for Epic naming, split on "_" instead?
        if ! [ $sample_name == $segment_sample_name ]; then
            echo "ERROR: sample names do not match between sample VCF and segment VCF"
            exit 1
        fi

        # identify CNVs from the segments VCF (bcftools filter)
        # parse their coordinates to a bed file (bcftools query)
        bcftools filter -e 'FORMAT/GT=="0/0"' $segments_vcf_path | \
            bcftools query -f '%CHROM\t%POS\t%INFO/END\n' > CNV_coords.bed

        # check whether any of the CNVs cover any PRS variants of interest
        bedtools intersect -a $PRS_variants_path -b CNV_coords.bed > CNV_intersect.bed

        # write list of affected PRS positions to output file
        echo "The following PRS variants are potentially found within a CNV. Please investigate further." > "$sample_name"_cnv_check.txt
        # with header from PRS bed
        grep "#" $PRS_variants_path >> "$sample_name"_cnv_check.txt
        cat CNV_intersect.bed >> "$sample_name"_cnv_check.txt
        echo -e "\nEnd of file" >> "$sample_name"_cnv_check.txt
    else
        echo "CNV checking was not performed for this sample." > "$sample_name"_cnv_check.txt
    fi

    ## 2. check coverage:
    mark-section "Checking for low coverage"
    # identify variants with low coverage (below depth threshold input)
    filter="FORMAT/DP<$depth"
    bcftools filter -i $filter $sample_vcf_path > low_cov_PRS.vcf

    # write list of affected PRS variants to output file
    echo "The following PRS variants are not covered to $depth x:" > "$sample_name"_coverage_check.txt
    echo -e "#CHROM\tPOS\tREF\tALT\tDP" >> "$sample_name"_coverage_check.txt
    # write relevant info about affected variants from filtered VCF
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\n' low_cov_PRS.vcf  >> "$sample_name"_coverage_check.txt
    echo -e "\nEnd of file" >> "$sample_name"_coverage_check.txt


    ## 3. convert VCF file:
    mark-section "Filtering PRS VCF"
    # exclude low covered PRS variants
    bcftools filter -e $filter $sample_vcf_path > "$sample_name"_canrisk_PRS.vcf
    grep -v ^# "$sample_name"_canrisk_PRS.vcf | wc -l

    ### OUTPUTS
    mark-section "Uploading output files"
    # upload files
    canrisk_VCF=$(dx upload "$sample_name"_canrisk_PRS.vcf --brief)
    cnv_check=$(dx upload "$sample_name"_cnv_check.txt --brief)
    coverage_check=$(dx upload "$sample_name"_coverage_check.txt --brief)

    dx-jobutil-add-output canrisk_PRS "$canrisk_VCF" --class=file
    dx-jobutil-add-output cnv_check "$cnv_check" --class=file
    dx-jobutil-add-output coverage_check "$coverage_check" --class=file

	dx-upload-all-outputs --parallel
	mark-success

}
