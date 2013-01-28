source ../setenv.sh
${mwa_reduce_dir}/calcgains -a 24 preprocessed.ms 207.701 -0.48 123.5 gains-avg.txt
${mwa_reduce_dir}/applygains preprocessed.ms gains-avg.txt
${mwa_reduce_dir}/calcgains preprocessed.ms 207.701 -0.48 123.5 gains2.txt
gnuplot gains2.plt
okular gains2.ps
