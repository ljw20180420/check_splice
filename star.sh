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

star_map_all() {
    local exp
    for exp in total pro clip
    do
        if [[ ${exp} == "total" ]]
        then
            local exp_dir="total_rna_seq"
        else
            local exp_dir="${exp}_seq"
        fi
        local protein
        for protein in $(ls ${root_dir}/${exp_dir})
        do
            if [[ ! -d "${root_dir}/${exp_dir}/${protein}" ]]
            then
                continue
            fi
            if [[ -d "${root_dir}/${exp_dir}/${protein}/non_rDNA" ]]
            then
                local data_dir="${root_dir}/${exp_dir}/${protein}/non_rDNA"
            else
                local data_dir="${root_dir}/${exp_dir}/${protein}"
            fi
            for R1 in $(ls ${data_dir}/*_R1.fq.gz)
            do
                local base=${R1##*/}
                local stem=${base%_R1.fq.gz}
                local clone rep
                IFS="_" read clone rep <<<${stem}

                if [[ -f "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam" ]]
                then
                    echo "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam exists"
                    continue
                fi

                echo "process ${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"

                star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                    "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                    "${data_dir}/${clone}_${rep}_R2.fq.gz"

                mv "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Aligned.sortedByCoord.out.bam" \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}.bam"
                rm -r \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}_STARgenome" \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.out" \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.progress.out" \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}Log.final.out" \
                    "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}SJ.out.tab"
            done
        done
    done
}

star_index="/home/ljw/.local/share/genomes/GRCh37/index/star"
root_dir="/home/ljw/sdc1/hush"
# index_ribosome
# filter_ribosome_all
# star_map_all
