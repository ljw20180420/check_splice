#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

check_existence() {
    local file=$1

    if [[ -s "${file}" ]]
    then
        echo "File exists"
        exit 1
    fi
}

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
    local protein_dir=$1
    local clone=$2
    local rep=$3

    bowtie2 --no-unal -I 1 -X 1000 -p 16 \
        -x rdnaModel/rdnaModel \
        -1 ${protein_dir}/${clone}_${rep}_R1.fq.gz \
        -2 ${protein_dir}/${clone}_${rep}_R2.fq.gz \
        --un-conc-gz ${protein_dir}/${clone}_${rep}.rDNAremoved.fq.gz \
        2> ${protein_dir}/${clone}_${rep}.rDNA.txt |
    samtools view -S -b -o ${protein_dir}/${clone}_${rep}.rDNA.bam -

    mkdir -p "${protein_dir}/non_rDNA"
    mv "${protein_dir}/${clone}_${rep}.rDNAremoved.fq.1.gz" \
        "${protein_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"
    mv "${protein_dir}/${clone}_${rep}.rDNAremoved.fq.2.gz" \
        "${protein_dir}/non_rDNA/${clone}_${rep}_R2.fq.gz"
    mv "${protein_dir}/${clone}_${rep}.rDNA.bam" \
        "${protein_dir}/${clone}_${rep}.rDNA.txt" \
        "${protein_dir}/non_rDNA/"
}

filter_ribosome_all() {
    local exp
    for exp in pro clip
    do
        local protein
        for protein in $(ls ${root_dir}/${exp}_seq)
        do
            if [[ ! -d "${root_dir}/${exp}_seq/${protein}" ]]
            then
                continue
            fi
            for R1 in $(ls ${root_dir}/${exp}_seq/${protein}/*_R1.fq.gz)
            do
                local base=${R1##*/}
                local stem=${base%_R1.fq.gz}
                local clone rep
                IFS="_" read clone rep <<<${stem}

                if [[ -s "${root_dir}/${exp}_seq/${protein}/non_rDNA/${clone}_${rep}_R1.fq.gz" ]]
                then
                    echo "${root_dir}/${exp}_seq/${protein}/non_rDNA/${clone}_${rep} exists"
                    continue
                fi
                echo "process ${root_dir}/${exp}_seq/${protein}/non_rDNA/${clone}_${rep}"

                filter_ribosome \
                    "${root_dir}/${exp}_seq/${protein}" \
                    "${clone}" \
                    "${rep}"
            done
        done
    done
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
        --outFilterIntronMotifs RemoveNoncanonical \
        --sjdbGTFfile ${root_dir}/data/hg19.ncbiRefSeq.gtf
}

clear_map() {
    local exp=$1
    local protein=$2
    local clone=$3
    local rep=$4
    
    mv "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Aligned.sortedByCoord.out.bam" \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"
    rm -r \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}_STARgenome" \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.out" \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.progress.out" \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.final.out" \
        "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}SJ.out.tab"
}

# star_map_all() {
#     for exp in total pro clip
#     do
#         if [[ ${exp} == "total" ]]
#         then
#     done
# }

total_map() {
    local exp="total"

    local protein="PPHLN1"
    local data_dir="${root_dir}/total_rna_seq/${protein}"

    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2" "rep3"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done
}

pro_map() {
    local exp="pro"

    local protein="MPP8"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D22" "DF15" "DF17" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D173" "DCF98" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="PPHLN1"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="TASOR"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D5" "DF4" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done
}

clip_map() {
    local exp="clip"

    local protein="MPP8"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "M4" "M5" "M15"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "C6"
    do
        for rep in "1-High1" "1-High2" "1-Low1" "1-Low2" "2-High1" "2-High2" "2-Low1" "2-Low2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "N77"
    do
        for rep in "1-Large-rep1" "1-Large-rep2" "1-Middle-rep1" "1-Middle-rep2" "1-Small-rep1" "1-Small-rep2" \
            "2-Large-rep1" "2-Large-rep2" "2-Middle-rep1" "2-Middle-rep2" "2-Small-rep1" "2-Small-rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="PPHLN1"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "PP98" "PP304"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="TASOR"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "TA242"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done

    local protein="WT"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "WT"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"

            clear_map ${exp} ${protein} ${clone} ${rep}
        done
    done
}

star_index="/home/ljw/.local/share/genomes/GRCh37/index/star"
root_dir="/home/ljw/sdc1/hush"
# index_ribosome
filter_ribosome_all
# total_map
# pro_map
# clip_map
