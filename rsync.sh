#!/bin/bash

source_dir="${HOME}/sdc1/hush"
target_dir="106:hush"
# the ending slash / matters
for folder in "result/hic/interact/" "bam/merge/precursor/" "bam/merge/splice/"
do
    rsync -avz --delete "${source_dir}/${folder}" "${target_dir}/${folder}"
done
