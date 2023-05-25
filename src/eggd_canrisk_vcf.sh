#!/bin/bash

set -exo pipefail

main() {
    ### INPUTS

    # get sample vcf
    # get PRS list
    # get intervals list for sample (if provided)

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    echo $PRS_variants_path $sample_vcf_path $intervals_file_path

    exit(0)

    ### PROCESS

    # make the bed file chr, start and end, then remove the first rows
    # which is the alpha value for the PRS risk and headers
    awk -F "," '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $PRS_variants_path | tail -n +3  > PRS_variants.bed

    # TODO - make the above generic. I.e. allow two inputs, a bed OR a PRS file. Convert PRS to bed here if needed, otherwise just use bed.

    # sort both files
    sort -k1V -k2n -k3n PRS_variants.bed > sorted_targets.bed 
    

    # if intervals provided
    #   check if any PRS positions are in a cnv (see process on other page)
    if $intervals_file; then
        # get rid of the header
        grep -v ^# intervals_file.vcf > intervals_no_header.vcf

        # make a list of any intervals that exhibit copy number variation
        rm found_cnvs.bed
        touch found_cnvs.bed
        while read line; do 
            chrom=$(echo $line | awk -F ' ' '{ print $3 }' | awk -F '_' '{ print $2 }');
            start=$(echo $line | awk -F ' ' '{ print $3 }' | awk -F '_' '{ print $3 }');
            end=$(echo $line | awk -F ' ' '{ print $3 }' | awk -F '_' '{ print $4 }');
            genotype=$(echo $line | awk -F ' ' '{ print $10 }' | awk -F ':' '{ print $1 }');
            if (($genotype > 0)); then
                echo -e $chrom"\t"$start"\t"$end >> found_cnvs.bed;
            fi;
        done < intervals_no_header.vcf
    
        # check whether any of the CNVs cover any PRS variants of interest
        bedtools intersect -a EGLH_303.bed -b found_cnvs.bed > intersect.bed

        # get proper details of any PRS variants found & alert the user via output file
        echo "The following PRS variants are found within a CNV. Please investigate further." > CNV_CHECK.txt
        while read line; do
        variant=$(echo $line | awk -F ' ' '{print $1","$3}');
        grep $variant $PRS_variants_path >> CNV_CHECK.txt;
        done < intersect.bed
    fi

    # modify genotypes to remove non-PRS alleles
    grep ^# sample_vcf > header_txt
    grep -v ^# sample_vcf > no_header

    while read line; do
        POSITION=$(echo $line | awk -F ' ' '{ print $1"\t"$2 }')
        PRS_ALT=$(grep $POSITION PRS_variants.bed | awk -F '\t' '{ print $5 }')
        if [[ echo $line | awk -F ' ' '{ print $10 }' | awk -F ':' '{ print $1 }' == '0/0']]; then
            output=$(echo line | awk -F ' ' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$PRS_ALT"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }')
    done < no_header

    # for each line (after header)
    #   if GT is 0/0:
    #     all is well. Do nothing/ add to output as is.
    #     not quite - need to add PRS alt allele
    #   else:
    #     if ALT == PRS_ALT:
    #       all is well. Do nothing/ add to output as is.
    #     elif ALT contains PRS_ALT (multi-allelic):
    #       extract PRS_ALT, dump the rest, set GT to unk/1

    ### OUTPUT

    # save as filename_canrisk.vcf
    # save cnv_report if needed
    # upload to project

    # upload VCF
    filtered_VCF=$(dx upload $out_file --brief)
    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file

}