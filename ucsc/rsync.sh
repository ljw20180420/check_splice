#!/bin/bash

# change to the dir of the script
cd $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

source_dir="${HOME}/sdc1/hush"
target_dir="106:hush"
# the ending slash / matters
for folder in "result/hic/interact/" "bam/merge/precursor/" "bam/merge/splice/"; do
	rsync -avz --delete "${source_dir}/${folder}" "${target_dir}/${folder}"
done
rsync -avz config.json 106:run_ucsc_genomebrowser_locally/src/config.json
