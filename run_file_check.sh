#!/usr/bin/env bash

for dir in $( ls -1 /mnt/)
do
    echo "Bioinformatics check for ${dir}"
    python3 file_check.py \
        --manifest /mnt/${dir}/manifest.txt \
        --datadir /mnt/${dir} \
        --outdir /home/output/${dir}_check \
        --logfile FILECHECK_LOG.log \
        --max_num_process 50 \
        --polling_period 1
done 2>/home/output/FILECHECK_STDERR_$(date '+%d%m%y_%H%M').log 1>/home/output/FILECHECK_STDOUT_$(date '+%d%m%y_%H%M').log
