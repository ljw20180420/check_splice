#!/bin/bash

draw_link() {
    local pairs=$1
    local region_key=$2
    local region=$3
    local output="${pairs%.pairs}.${region_key}.pdf"

    awk -F $'\t' -v OFS=$'\t' '
        NR > 1 {
            print $2, $3 - 1, $3, $4, $5 - 2, $5 - 1
        }
    ' ${pairs} |
    csvtk freq -tHf1,2,3,4,5,6 \
        > temp.links

    make_tracks_file \
        --trackFiles \
            temp.links \
            cpcdh.bed \
            pCBS.bed \
        -o tracks.ini

    sed -E -i '
        0,/^color = red$/{
            s/^color = red$/color = Reds/
        }
    ' tracks.ini

    sed -E -i '
        0,/^#min_value = 0$/{
            s/^#min_value = 0$/min_value = 0/
        }
    ' tracks.ini

    sed -E -i '
        0,/^#max_value = 1.2$/{
            s/^#max_value = 1.2$/max_value = '${max_value}'/
        }
    ' tracks.ini

    sed -E -i '
        0,/^#line_width = 0.5$/{
            s/^#line_width = 0.5$/line_width = 0.5/
        }
    ' tracks.ini

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
max_value=50.0


csvtk filter2 -f '$type == "exon"' "${root_dir}/result/cpcdh.csv" |
csvtk cut -f chrom,start,end,name |
csvtk del-header -T \
    > cpcdh.bed

csvtk cut -tf1-4 \
    "${root_dir}/result/pCBS.bed" \
    > pCBS.bed

for key in cpcdh alpha beta gamma
do
    region="${regions[${key}]}"
    for orientation in ff rr
    do
        draw_link \
            ${root_dir}/result/hic/${orientation}.pairs \
            ${key} ${region} 2> error.log
        exit 0
        while read exp_protein_wt
        do
            draw_link \
                ${root_dir}/result/hic/${exp_protein_wt}_${orientation}.pairs \
                ${key} ${region}
        done < exp_protein_wts.txt
    done
done
