#!/bin/bash

set -exo pipefail

main() {
    ### INPUTS

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    ### PROCESS

    # make the bed file chr, start and end, then remove the first rows
    # which is the alpha value for the PRS risk and headers
    awk -F "," '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $PRS_variants_path | tail -n +3  > PRS_variants.bed

    # TODO - make the above generic. I.e. allow two inputs, a bed OR a PRS file. Convert PRS to bed here if needed, otherwise just use bed.

    # if intervals provided
    #   check if any PRS positions are in a cnv (see process on other page)
    if [ $intervals_file_path ]; then
        # get rid of the header
        grep -v ^# $intervals_file_path > intervals_no_header.vcf
        # make a list of any intervals that exhibit copy number variation
        touch found_cnvs.bed
        while read line; do 
            chrom=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $2 }');
            start=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $3 }');
            end=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $4 }');
            genotype=$(printf '%s\t' $line | awk -F '\t' '{ print $10 }' | awk -F ':' '{ print $1 }');
            # mprint relevant lines to file (GT = 1 or 2 if DEL/DUP)
            if [ $genotype > 0 ]; then
                printf '%s\t%s\t%s\n' $chrom $start $end >> found_cnvs.bed;
            fi;
        done < intervals_no_header.vcf
        # check whether any of the CNVs cover any PRS variants of interest
        bedtools intersect -a PRS_variants.bed -b found_cnvs.bed > intersect.bed
        # get proper details of any PRS variants found & alert the user via output file
        echo "The following PRS variants are found within a CNV. Please investigate further." > CNV_CHECK.txt
        while read line; do
            variant=$(printf '%s\t' $line | awk -F '\t' '{print $1","$3}');
            grep $variant $PRS_variants_path >> CNV_CHECK.txt;
        done < intersect.bed
    fi

    # norm/decompose vcf to split multi-allelics
    bcftools norm -m- $sample_vcf_path > vcf_norm

    # grab header
    grep ^# vcf_norm > canrisk_VCF
    grep -v ^# vcf_norm > no_header

    # make modified vcf
    while read line; do
        echo $line
        POSITION=$(printf '%s\t' $line | awk -F '\t' '{ print $1"\t"$3 }')
        REF=$(printf '%s\t' $line | awk -F '\t' '{ print $4 }')
        ALT=$(printf '%s\t' $line | awk -F '\t' '{ print $5 }')
        echo $POSITION $REF $ALT
        # TODO - add check for whether variant is in sample? or maybe just grep for position as ref might not be the same for 0/0 indels
        SAMPLE_LINE=$(grep "$POSITION$(printf '\t').*$(printf '\t')$REF" vcf_norm)
        GENOTYPE=$(printf '%s\t' $SAMPLE_LINE | awk -F '\t' '{ print $10 }' | awk -F ':' '{ print $1 }')
        echo $POSITION $REF $ALT $SAMPLE_LINE $GENOTYPE
        # check if grep finds the variant
        if grep -q "$POSITION$(printf '\t').*$(printf '\t')$REF$(printf '\t')$ALT" vcf_norm; then
            var = $(grep -q "$POSITION$(printf '\t').*$(printf '\t')$REF$(printf '\t')$ALT" vcf_norm)
            output=$(printf '%s\t' $var)
        # if not then check if genotype is 0/0 and add in PRS ALT
        elif $GENOTYPE == '0/0'; then
            echo $GENOTYPE
            output = $(printf '%s\t' $SAMPLE_LINE | awk -v var="$ALT" -F '\t' '{ print $1"\t"$2"\t"$3"\t"$4"\t"var"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }')
        # otherwise the VCF is wrong as we should be force genotyping all PRS positions
        else
            echo PRS variant not found in sample VCF. Please check it has been generated correctly.
            continue
        fi
        # print relevant lines to output file
        echo $output
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $output >> canrisk_VCF
    done < PRS_variants.bed

    # upload VCF
    canrisk_VCF=$(dx upload canrisk_VCF --brief)
    cnv_check=$(dx upload CNV_CHECK.txt --brief)
    dx-jobutil-add-output canrisk_VCF "$canrisk_VCF" --class=file
    dx-jobutil-add-output CNV_CHECK "$cnv_check" --class=file

}