#!/bin/bash

set -exo pipefail

main() {
    ### INPUTS

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    ### PREPARE
    mark-section "Preparing reference genome files"
    # move reference genome files to /home/dnanexus/
    mv $reference_fasta_path /home/dnanexus/ref_gen.fa.gz
    mv $reference_fasta_index_path /home/dnanexus/ref_gen.fa.fai
    # unpack reference genome fasta
    gunzip /home/dnanexus/ref_gen.fa.gz

    # # zip and fasta reference properly
    # zcat $reference_fasta_path | bgzip -c > ref.bgz

    # normalise/decompose VCF to split multi-allelics
    mark-section "Normalising input VCF"
    bcftools norm -m- -f ref_gen.fa -cs $sample_vcf_path > norm.vcf

    # get sample name
    sample_name=$(basename $sample_vcf_path | awk -F '-' '{ print $1"-"$2 }')
    # TODO future proofing for Epic naming, split on "_" instead?


    ### PROCESS

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
    # identify variants with low coverage (<20 DP)
    # parse their coordinates to a bed file (bcftools query)
    bcftools filter -i 'INFO/DP<20' norm.vcf | \
        bcftools query -f '%CHROM\t%POS\t%INFO/END\n' > low_cov.bed

    # check whether any of the low covered positions overlap any PRS variants of interest
    bedtools intersect -a $PRS_variants_path -b low_cov.bed > cov_intersect.bed

    # write list of affected PRS positions to output file
    echo "The following PRS positions are not covered to 20x:" > "$sample_name"_coverage_check.txt
    # with header from PRS bed
    grep "#" $PRS_variants_path >> "$sample_name"_coverage_check.txt
    cat cov_intersect.bed >> "$sample_name"_coverage_check.txt
    # TODO need some logic to check whether there are any
    # and if so, add DP to each variant from VCF
    # echo -e "# CHROM\tPOS\tDP" >> "$sample_name"_coverage_check.txt
    echo -e "\nEnd of file" >> "$sample_name"_coverage_check.txt


    ## 3. convert VCF file:
    mark-section "Converting input VCF"
    # exclude low covered PRS variants?
    bcftools filter -e 'INFO/DP<20' norm.vcf > norm_cov.vcf

    # grab header
    grep ^# norm_cov.vcf > "$sample_name"_canrisk_PRS.vcf

    grep -v ^# $PRS_variants_path > PRS.bed

    # make modified vcf
    while read line; do
        POSITION=$(printf '%s\t' $line | awk -F '\t' '{ print $1"\t"$3 }')
        REF=$(printf '%s\t' $line | awk -F '\t' '{ print $4 }')
        ALT=$(printf '%s\t' $line | awk -F '\t' '{ print $5 }')
        # check position actually exists & warn if not (all positions should be present)
        if ! grep -qP "^$POSITION\t" norm_cov.vcf; then
            echo "ERROR: Position $POSITION not found in sample vcf - please ensure it has been formed correctly."
            exit 1
        fi
        SAMPLE_LINE=$(grep -P "^$POSITION\t" norm_cov.vcf)
        # get GT & DP indices
        FORMAT=$(printf '%s\t' $SAMPLE_LINE | awk -F '\t' '{ print $9 }')
        IFS=':' read -ra FIELD <<< "$FORMAT"
        idx=1
        for i in "${FIELD[@]}"; do
            if [ $i == "GT" ]; then GT_INDEX=$idx
            elif [ $i == "DP" ]; then DP_INDEX=$idx
            fi
            let "idx = idx + 1"
        done
        # get GT & DP values
        FORMAT_VALUES=$(printf '%s\t' $SAMPLE_LINE | awk -F '\t' '{ print $10 }')
        GENOTYPE=$(awk -v var="$GT_INDEX" -vt="${FORMAT_VALUES[*]}" 'BEGIN{n=split(t,a,":"); print a[var]}')
        DEPTH=$(awk -v var="$DP_INDEX" -vt="${FORMAT_VALUES[*]}" 'BEGIN{n=split(t,a,":"); print a[var]}')
        # # check depth is over 20x & record if not
        # if [ $DEPTH -lt 20 ]; then
        #     echo -e $POSITION "\t" $DEPTH >> "$sample_name"_coverage_check.txt
        #     # TODO output all info about PRS from bed, include header
        #     # TODO DITCH/remove/exclude records where depth < 20x
        # fi
        # check if grep finds the variant
        if grep -qP "^$POSITION\t.*\t$REF\t$ALT\t" norm_cov.vcf; then
            var=$(grep -P "^$POSITION\t.*\t$REF\t$ALT\t" norm_cov.vcf)
            output=$(printf '%s\t' $var)
            echo "variant is in norm_cov"
            echo $output
        # if not then check if genotype is 0/0 and add in PRS ALT
        elif [ $GENOTYPE == '0/0' ]; then
            output=$(printf '%s\t' $SAMPLE_LINE | awk -v var1="$REF" -v var2="$ALT" -F '\t' '{ print $1"\t"$2"\t"$3"\t"var1"\t"var2"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }')
            echo "variant is NOT in norm_cov, and has GT==0/0"
            echo $output
        # otherwise PRS is not present so we recode as 0/0
        else
            output=$(printf '%s\t' $SAMPLE_LINE | awk -v var1="$REF" -v var2="$ALT" -v var3="$DEPTH" -F '\t' '{ print $1"\t"$2"\t.\t"var1"\t"var2"\t.\t.\t.\tGT:DP\t0/0:"var3 }')
            echo "variant is NOT in norm_cov, and NOT GT==0/0"
            echo "PRS variant not found in sample VCF. Adding line to say sample is 0/0 at this position."
            echo $output
        fi
        # print relevant lines to output file
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $output >> "$sample_name"_canrisk_PRS.vcf
    done < PRS.bed

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
