source ../setenv.sh
${mwa_reduce_dir}/calcphaseoffsets -a 32 -f preprocessed.ms phases.txt
gnuplot phases.plt
okular phases.ps
