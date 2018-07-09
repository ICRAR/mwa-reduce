#! /bin/bash
source ../setenv.sh
${mwa_reduce_dir}/flagmwa preprocessed.ms 0 24 1 4
${mwa_reduce_dir}/flagantennae preprocessed.ms 0 6 8 16 31
#${mwa_reduce_dir}/flagsubbands preprocessed.ms 24 0
