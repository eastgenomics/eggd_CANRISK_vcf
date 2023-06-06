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

    # get sample name
    sample_name=$(basename $sample_vcf_path | awk -F '-' '{ print $1"-"$2 }')

    # if segments provided
    #   check if any PRS positions are in a cnv (see process on other page)
    if [ $segments_file_path ]; then
        # check sample matches vcf
        segment_sample_name=$(basename $segments_file_path | awk -F '-' '{ print $1"-"$2 }')
        if ! [ $sample_name == $segment_sample_name ]; then
            echo "ERROR: sample names do not match between sample vcf and segment vcf"
            exit 1
        fi
        # get rid of the header
        grep -v ^# $segments_file_path > segments_no_header.vcf
        # make a list of any intervals that exhibit copy number variation
        touch found_cnvs.bed
        while read line; do 
            chrom=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $2 }');
            start=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $3 }');
            end=$(printf '%s\t' $line | awk -F '\t' '{ print $3 }' | awk -F '_' '{ print $4 }');
            genotype=$(printf '%s\t' $line | awk -F '\t' '{ print $10 }' | awk -F ':' '{ print $1 }');
            echo $chrome $start $end $genotype
            # print relevant lines to file (GT = 1 or 2 if DEL/DUP)
            if ! [ $genotype == '0/0' ]; then
                printf '%s\t%s\t%s\n' $chrom $start $end >> found_cnvs.bed;
            fi;
        done < segments_no_header.vcf
        # check whether any of the CNVs cover any PRS variants of interest
        bedtools intersect -a PRS_variants.bed -b found_cnvs.bed > intersect.bed
        # get proper details of any PRS variants found & alert the user via output file
        echo "The following PRS variants are potentially found within a CNV. Please investigate further." > cnv_check.txt
        while read line; do
            variant=$(printf '%s\t' $line | awk -F '\t' '{print $1","$3}');
            grep $variant $PRS_variants_path >> cnv_check.txt;
        done < intersect.bed
        echo -e "\nEnd of file" >> cnv_check.txt
    else
        echo "CNV checking was not requested for this sample." > cnv_check.txt
    fi

    # norm/decompose vcf to split multi-allelics
    bcftools norm -m- $sample_vcf_path > vcf_norm

    # grab header
    grep ^# vcf_norm > canrisk_VCF

    # initiate coverage file
    echo "The following PRS positions are not covered to 20x:" > coverage_check.txt
    echo -e "# CHROM\tPOS\tDP" >> coverage_check.txt

    # make modified vcf
    while read line; do
        POSITION=$(printf '%s\t' $line | awk -F '\t' '{ print $1"\t"$3 }')
        REF=$(printf '%s\t' $line | awk -F '\t' '{ print $4 }')
        ALT=$(printf '%s\t' $line | awk -F '\t' '{ print $5 }')
        # check position actually exists & warn if not (all positions should be present)
        if ! [ grep "^$POSITION$(printf '\t')" vcf_norm ]; then
            echo "ERROR: Position $POSITION not found in sample vcf - please ensure it has been formed correctly."
            exit 1
        fi
        SAMPLE_LINE=$(grep "^$POSITION$(printf '\t')" vcf_norm)
        GENOTYPE=$(printf '%s\t' $SAMPLE_LINE | awk -F '\t' '{ print $10 }' | awk -F ':' '{ print $1 }')
        DEPTH=$(printf '%s\t' $SAMPLE_LINE | awk -F '\t' '{ print $10 }' | awk -F ':' '{ print $3 }')
        # check depth is over 20x & record if not
        if [ $DEPTH -lt 20 ]; then
            echo -e $POSITION "\t" $DEPTH >> coverage_check.txt
        fi
        # check if grep finds the variant
        if grep -q "^$POSITION$(printf '\t').*$(printf '\t')$REF$(printf '\t')$ALT$(printf '\t')" vcf_norm; then
            var=$(grep "^$POSITION$(printf '\t').*$(printf '\t')$REF$(printf '\t')$ALT$(printf '\t')" vcf_norm)
            output=$(printf '%s\t' $var)
        # if not then check if genotype is 0/0 and add in PRS ALT
        elif [ $GENOTYPE == '0/0' ]; then
            output=$(printf '%s\t' $SAMPLE_LINE | awk -v var="$ALT" -F '\t' '{ print $1"\t"$2"\t"$3"\t"$4"\t"var"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }')
        # otherwise the VCF is wrong as we should be force genotyping all PRS positions
        else
            output=$(printf '%s\t' $SAMPLE_LINE | awk -v var1="$REF" -v var2="$ALT" -v var3="$DEPTH" -F '\t' '{ print $1"\t"$2"\t.\t"var1"\t"var2"\t.\t.\t.\tGT:DP\t0/0:"var3 }')
            echo "PRS variant not found in sample VCF. Adding line to say sample is 0/0 at this position."
        fi
        # print relevant lines to output file
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $output >> canrisk_PRS.vcf
    done < PRS_variants.bed
    echo -e "\nEnd of file" >> coverage_check.txt

    # upload VCF
    canrisk_VCF=$(dx upload canrisk_PRS.vcf --brief)
    cnv_check=$(dx upload cnv_check.txt --brief)
    coverage_check=$(dx upload coverage_check.txt --brief)
    dx-jobutil-add-output canrisk_PRS "$canrisk_VCF" --class=file
    dx-jobutil-add-output cnv_check "$cnv_check" --class=file
    dx-jobutil-add-output coverage_check "$coverage_check" --class=file

}
