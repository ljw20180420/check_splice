#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

index_ribosome() {
    mkdir -p rdnaModel
    bigBedToBed https://hgdownload.soe.ucsc.edu/gbdb/hs1/rdnaModel/rdnaModel.bb rdnaModel/rdnaModel.bed
    bedtools getfasta \
        -fi /home/ljw/.local/share/genomes/hs1/hs1.fa \
        -bed rdnaModel/rdnaModel.bed \
        -fo rdnaModel/rdnaModel.fa
    bowtie2-build rdnaModel/rdnaModel.fa rdnaModel/rdnaModel
}

filter_ribosome() {
    local R1=$1
    local R2=$2
    local output=$3

    bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p 16 \
        -x rdnaModel/rdnaModel \
        -1 ${R1} -2 ${R2} \
        --un-conc-gz ${output}.rDNAremoved.fq.gz \
        2> ${output}.rDNA.txt |
    samtools view -S -b -o ${output}.rDNA.bam - 
}

pro_filter_ribosome() {
    local data_dir="${root_dir}/pro_seq/MPP8"
    mkdir -p ${data_dir}/non_rDNA
    for clone in "D22" "DF15" "DF17" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            if [[ -s "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz" ]]
            then
                continue
            fi

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            mv "${data_dir}/${clone}_${rep}.rDNAremoved.fq.1.gz" \
                ${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz
            mv "${data_dir}/${clone}_${rep}.rDNAremoved.fq.2.gz" \
                ${data_dir}/non_rDNA/${clone}_${rep}_R2.fq.gz
            mv "${data_dir}/${clone}_${rep}.rDNA.bam" \
                "${data_dir}/${clone}_${rep}.rDNA.txt" \
                "${data_dir}/non_rDNA/"
        done
    done

    local data_dir="${root_dir}/pro_seq/NP229"
    for clone in "D173" "DCF98" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"
        done
    done

    local data_dir="${root_dir}/pro_seq/PPHLN1"
    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"
        done
    done

    local data_dir="${root_dir}/pro_seq/TASOR"
    for clone in "D5" "DF4" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"
        done
    done
}

get_common_prefix() {
    # If no arguments provided, exit
    [ $# -eq 0 ] && return

    local prefix="$1"
    shift

    # Compare prefix against all other arguments
    for item in "$@"; do
        # Truncate prefix until it matches the start of the item
        while [[ "$item" != "$prefix"* ]]; do
            prefix="${prefix%?}"
            [ -z "$prefix" ] && break 2
        done
    done

    echo "$prefix"
}

star_map() {
    local prefix=$1
    shift

    STAR \
        --runThreadN 16 \
        --genomeDir ${star_index} \
        --readFilesIn "$@" \
        --readFilesCommand zcat \
        --outFileNamePrefix ${prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --sjdbGTFfile ${root_dir}/hg19.ncbiRefSeq.gtf
}

total_PPHLN1_map() {
    for clone in "PPD169" "PPD221" "WT6"
    do
        for rep in "rep1" "rep2" "rep3"
        do
            star_map "${root_dir}/bam/total_PPHLN1_${clone}_${rep}" \
                "${root_dir}/total_rna_seq/PPHLN1/${clone}_Total_${rep}.R1.raw.fastq.gz" \
                "${root_dir}/total_rna_seq/PPHLN1/${clone}_Total_${rep}.R2.raw.fastq.gz"
        done
    done
}

pro_TASOR_map() {
    for clone in "TA-D5" "TA-DF4" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/pro_TASOR_${clone}_${rep}" \
                "${root_dir}/pro_seq/TASOR/${clone}_${rep}_R1.fastq.gz" \
                "${root_dir}/pro_seq/TASOR/${clone}_${rep}_R2.fastq.gz"
        done
    done
}

pro_MPP8_map() {
    for clone in "M8-D22" "M8-DF15" "M8-DF17" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/pro_MPP8_${clone}_${rep}" \
                "${root_dir}/pro_seq/MPP8/${clone}_${rep}_R1.fastq.gz" \
                "${root_dir}/pro_seq/MPP8/${clone}_${rep}_R2.fastq.gz"
        done
    done
}

pro_NP220_map() {
    for clone in "NP220-D173" "NP220-DCF98" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/pro_NP220_${clone}_${rep}" \
                "${root_dir}/pro_seq/NP220/${clone}_${rep}_R1.fastq.gz" \
                "${root_dir}/pro_seq/NP220/${clone}_${rep}_R2.fastq.gz"
        done
    done
}

clip_WT_map() {
    for clone in "WT"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/clip_WT_${clone}_${rep}" \
                "${root_dir}/clip_seq/WT/${clone}_${rep}_R1.fq.gz" \
                "${root_dir}/clip_seq/WT/${clone}_${rep}_R2.fq.gz"
        done
    done
}

clip_MPP8_map() {
    for clone in "M4" "M5" "M15"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/clip_MPP8_${clone}_${rep}" \
                "${root_dir}/clip_seq/MPP8/${clone}_${rep}_R1.fq.gz" \
                "${root_dir}/clip_seq/MPP8/${clone}_${rep}_R2.fq.gz"
        done
    done
}

clip_TASOR_map() {
    for clone in "TA242"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/clip_TASOR_${clone}_${rep}" \
                "${root_dir}/clip_seq/TASOR/${clone}_${rep}.R1.raw.fastq.gz" \
                "${root_dir}/clip_seq/TASOR/${clone}_${rep}.R2.raw.fastq.gz"
        done
    done
}

clip_PPHLN1_map() {
    for clone in "PP98" "PP304"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/clip_PPHLN1_${clone}_${rep}" \
                "${root_dir}/clip_seq/PPHLN1/${clone}_${rep}.R1.raw.fastq.gz" \
                "${root_dir}/clip_seq/PPHLN1/${clone}_${rep}.R2.raw.fastq.gz"
        done
    done
}

star_index="/home/ljw/.local/share/genomes/GRCh37/index/star"
root_dir="/home/ljw/sdc1/hush"
pro_filter_ribosome
