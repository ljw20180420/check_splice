#!/bin/bash

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
    local star_index=$1
    local prefix=$2
    shift 2

    STAR \
        --runThreadN 16 \
        --genomeDir ${star_index} \
        --readFilesIn "$@" \
        --readFilesCommand zcat \
        --outFileNamePrefix $(get_common_prefix TASORdel_Total_R1.fastq.gz TASORdel_Total_R2.fastq.gz) \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outWigType bedGraph \
        --outFilterIntronMotifs RemoveNoncanonical \
        --quantMode GeneCounts   
}

star_index=/home/ljw/.local/share/genomes/GRCh37/index/star
star_map "${star_index}" "$@"
