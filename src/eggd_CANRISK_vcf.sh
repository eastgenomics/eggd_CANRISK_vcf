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
    # if segments file input is provided, check if any PRS variants are in a CNV
    # by converting CNV coordinates to a bed file
    if [ $segments_vcf_path ]; then
        mark-section "Checking for overlapping CNVs"
        # check CNV VCF filename matches sample from input VCF
        segment_sample_name=$(basename $segments_vcf_path | awk -F '_' '{ print $1 }')
        if ! [ $sample_name == $segment_sample_name ]; then
            echo "ERROR: sample names do not match between sample VCF and segment VCF"
            exit 1
        fi

        # identify CNVs from the segments VCF (bcftools filter)
        # parse their coordinates to a bed file (bcftools query)
        bcftools filter -e 'FORMAT/GT=="0/0"' $segments_vcf_path | \
            bcftools query -f '%CHROM\t%POS\t%INFO/END\n' > CNV_coords.bed

        # check whether any of the CNVs cover any PRS variants of interest
        bedtools intersect -a $sample_vcf_path -b CNV_coords.bed > PRS_intersect_CNV.tsv

        # write list of affected PRS variants to output file
        echo "The following PRS variants are located within a CNV call in this sample:" > "$sample_name"_cnv_check.txt
        echo -e "\n#CHROM\tPOS\tREF\tALT" >> "$sample_name"_cnv_check.txt
        # write relevant info about affected variants
        cut -f 1,2,4,5 PRS_intersect_CNV.tsv >> "$sample_name"_cnv_check.txt
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
    echo "The following PRS variants are not covered to $depth x read depth:" > "$sample_name"_coverage_check.txt
    echo -e "\n#CHROM\tPOS\tREF\tALT\tDP" >> "$sample_name"_coverage_check.txt
    # write relevant info about affected variants from filtered VCF
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\n' low_cov_PRS.vcf  >> "$sample_name"_coverage_check.txt
    echo -e "\nEnd of file" >> "$sample_name"_coverage_check.txt

    ## 3. filter VCF file:
    # TODO make exclusion optional
    mark-section "Filtering PRS VCF"
    # exclude low covered PRS variants
    bcftools filter -e $filter $sample_vcf_path > "$sample_name"_canrisk_PRS.vcf
    echo "Retained $(grep -v ^# "$sample_name"_canrisk_PRS.vcf | wc -l) variants after depth filtering"

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
