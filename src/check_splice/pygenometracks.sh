#!/bin/bash

draw_hic() {
    local hic=$1
    local region_key=$2
    local region=$3
    local output="${hic%.hic}.${region_key}.pdf"

    make_tracks_file \
        --trackFiles \
            ${hic} \
            <(
                tail -n+2 "${root_dir}/result/cpcdh.csv" |
                cut -d, -f1-6 |
                tr ',' '\t'
            ) \
            "${root_dir}/result/pCBS.bed" \
        -o tracks.ini
    pyGenomeTracks \
        --tracks tracks.ini \
        --region ${region} \
        --outFileName ${output}
}

root_dir="/home/ljw/sdc1/hush"
declare -A regions
regions[cpcdh]="chr5:140,125,000-140,922,000"
regions[alpha]="chr5:140,125,000-140,414,000"
regions[beta]="chr5:140,414,000-140,668,000"
regions[gamma]="chr5:140,668,000-140,922,000"

for key in cpcdh alpha beta gamma
do
    region="${regions[${key}]}"
    for orientation in ff rr
    do
        draw_hic \
            ${root_dir}/result/hic/${orientation}.hic \
            ${key} ${region}
        while read exp_protein_wt
        do
            draw_hic \
                ${root_dir}/result/hic/${exp_protein_wt}_${orientation}.hic \
                ${key} ${region}
        done < src/check_splice/exp_protein_wts.txt
    done
done
