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
    local R1=$1
    local R2=$2
    local output=$3

    bowtie2 --no-unal -I 1 -X 1000 -p 16 \
        -x rdnaModel/rdnaModel \
        -1 ${R1} -2 ${R2} \
        --un-conc-gz ${output}.rDNAremoved.fq.gz \
        2> ${output}.rDNA.txt |
    samtools view -S -b -o ${output}.rDNA.bam - 
}

move_filter_ribosome() {
    local data_dir=$1
    local clone=$2
    local rep=$3

    mkdir -p ${data_dir}/non_rDNA
    mv "${data_dir}/${clone}_${rep}.rDNAremoved.fq.1.gz" \
        ${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz
    mv "${data_dir}/${clone}_${rep}.rDNAremoved.fq.2.gz" \
        ${data_dir}/non_rDNA/${clone}_${rep}_R2.fq.gz
    mv "${data_dir}/${clone}_${rep}.rDNA.bam" \
        "${data_dir}/${clone}_${rep}.rDNA.txt" \
        "${data_dir}/non_rDNA/"
}

pro_filter_ribosome() {
    local data_dir="${root_dir}/pro_seq/MPP8"
    for clone in "D22" "DF15" "DF17" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/pro_seq/NP220"
    for clone in "D173" "DCF98" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/pro_seq/PPHLN1"
    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/pro_seq/TASOR"
    for clone in "D5" "DF4" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done
}

clip_filter_ribosome() {
    local data_dir="${root_dir}/clip_seq/MPP8"
    for clone in "M4" "M5" "M15"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/clip_seq/NP220"
    for clone in "C6"
    do
        for rep in "1-High1" "1-High2" "1-Low1" "1-Low2" "2-High1" "2-High2" "2-Low1" "2-Low2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done
    for clone in "N77"
    do
        for rep in "1-Large-rep1" "1-Large-rep2" "1-Middle-rep1" "1-Middle-rep2" "1-Small-rep1" "1-Small-rep2" \
            "2-Large-rep1" "2-Large-rep2" "2-Middle-rep1" "2-Middle-rep2" "2-Small-rep1" "2-Small-rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/clip_seq/PPHLN1"
    for clone in "PP98" "PP304"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/clip_seq/TASOR"
    for clone in "TA242"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
        done
    done

    local data_dir="${root_dir}/clip_seq/WT"
    for clone in "WT"
    do
        for rep in "rep1" "rep2"
        do
            check_existence "${data_dir}/non_rDNA/${clone}_${rep}_R1.fq.gz"

            filter_ribosome \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz" \
                "${data_dir}/${clone}_${rep}"

            move_filter_ribosome ${data_dir} ${clone} ${rep}
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
        --sjdbGTFfile ${root_dir}/hg19.ncbiRefSeq.gtf
}

total_map() {
    local exp="total"

    local protein="PPHLN1"
    local data_dir="${root_dir}/total_rna_seq/${protein}"

    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2" "rep3"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
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
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D173" "DCF98" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="PPHLN1"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D169" "D221" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="TASOR"
    local data_dir="${root_dir}/pro_seq/${protein}/non_rDNA"
    for clone in "D5" "DF4" "WT6"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
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
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "C6"
    do
        for rep in "1-High1" "1-High2" "1-Low1" "1-Low2" "2-High1" "2-High2" "2-Low1" "2-Low2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="NP220"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "N77"
    do
        for rep in "1-Large-rep1" "1-Large-rep2" "1-Middle-rep1" "1-Middle-rep2" "1-Small-rep1" "1-Small-rep2" \
            "2-Large-rep1" "2-Large-rep2" "2-Middle-rep1" "2-Middle-rep2" "2-Small-rep1" "2-Small-rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="PPHLN1"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "PP98" "PP304"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="TASOR"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "TA242"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done

    local protein="WT"
    local data_dir="${root_dir}/clip_seq/${protein}/non_rDNA"
    for clone in "WT"
    do
        for rep in "rep1" "rep2"
        do
            star_map "${root_dir}/bam/${exp}_${protein}_${clone}_${rep}" \
                "${data_dir}/${clone}_${rep}_R1.fq.gz" \
                "${data_dir}/${clone}_${rep}_R2.fq.gz"
        done
    done
}

star_index="/home/ljw/.local/share/genomes/GRCh37/index/star"
root_dir="/home/ljw/sdc1/hush"
# pro_filter_ribosome
# clip_filter_ribosome
total_map
pro_map
clip_map
