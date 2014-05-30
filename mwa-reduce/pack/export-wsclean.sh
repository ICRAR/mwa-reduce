#! /bin/bash
if [[ $1 == "" ]] ; then
		echo Syntax: $0 \<dir\>
else
		WORKDIR=`pwd`
		DEST="$1"
		cd ..
		cp -v areaset.* banddata.* beamevaluator.* buffered_lane.* cachedimageset.* fftconvolver.* fftresampler.* fitsreader.* fitswriter.* imagecoordinates.* imagebufferallocator.* imageweights.* inversionalgorithm.* lane.* layeredimager.* multibanddata.* nlplfitter.* matrix2x2.* model.* modelrenderer.* modelsource.* msselection.* polarizationenum.* progressbar.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* stopwatch.* uvector.* uvwdistribution.* weightmode.* wsclean.* wscleanmain.cpp wsinversion.* ${DEST}/
		cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h aocommon/uvector_03.h ${DEST}/aocommon/
		cp -v beam/lnaimpedance.* beam/tilebeam*.{h,cpp} beam/tileimpedance.{h,cpp} ${DEST}/beam/
		cp -v cleanalgorithms/*.{h,cpp} ${DEST}/cleanalgorithms/
		cp -v msprovider/*.{h,cpp} ${DEST}/msprovider/
		cp -v parser/*.h ${DEST}/parser/		
fi
