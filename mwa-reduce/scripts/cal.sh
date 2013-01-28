source ../setenv.sh
${mwa_reduce_dir}/phasecal -a 32 model-hyda-fieldE.txt preprocessed.ms phases.txt gains.txt
#${mwa_reduce_dir}/phasecal -a 32 hyda-141-model.txt preprocessed.ms phases.txt gains.txt
gnuplot phases.plt
gnuplot gains.plt
okular phases.ps gains.ps &

