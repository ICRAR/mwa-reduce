source ../setenv.sh
rm -f model-clean.txt model-clean-hires.txt
for (( a=0 ; a!=120 ; a++ )) ; do
    ${mwa_reduce_dir}/aoimage -scale 0.005 128 clean.fits preprocessed.ms
    cp clean.fits images/clean-iter-$a.fits

    echo Searching peak...
    ${mwa_reduce_dir}/imgclean clean.fits | tee model-clcomponent.txt

    echo Measuring spectrum of peak...
    ${mwa_reduce_dir}/spectrum model-clcomponent.txt preprocessed.ms > model-clcomp-spectrum.txt
    cat model-clcomp-spectrum.txt >> model-clean-hires.txt
    rm model-clcomponent.txt

    echo Interpolating spectrum...
    ${mwa_reduce_dir}/sdf -s 0.1 24 model-clcomp-spectrum.txt preprocessed.ms > model-clcomp-sdf.txt
    rm model-clcomp-spectrum.txt

    cat model-clcomp-sdf.txt >> model-clean.txt

    echo Subtracting component...
    ${mwa_reduce_dir}/subtrmodel model-clcomp-sdf.txt preprocessed.ms
    rm model-clcomp-sdf.txt
done
