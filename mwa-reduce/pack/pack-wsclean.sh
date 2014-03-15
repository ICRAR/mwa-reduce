WORKDIR=`pwd`
rm -rf /tmp/wsclean
mkdir /tmp/wsclean
mkdir /tmp/wsclean/aocommon
mkdir /tmp/wsclean/cleanalgorithms
mkdir /tmp/wsclean/msprovider
mkdir /tmp/wsclean/parser
cp CMakeLists-wsclean.txt /tmp/wsclean/CMakeLists.txt
cd ..
cp -v areaset.* banddata.* beamevaluator.* buffered_lane.* cachedimageset.* fftresampler.* fitsreader.* fitswriter.* imagecoordinates.* imagebufferallocator.* imageweights.* inversionalgorithm.* lane.* layeredimager.* multibanddata.* nlplfitter.* matrix2x2.* model.* modelrenderer.* modelsource.* msselection.* polarizationenum.* progressbar.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* stopwatch.* tilebeam.* uvector.* uvwdistribution.* weightmode.* wsclean.* wscleanmain.cpp wsinversion.* /tmp/wsclean/
cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h /tmp/wsclean/aocommon/
cp -v cleanalgorithms/*.{h,cpp} /tmp/wsclean/cleanalgorithms/
cp -v msprovider/*.{h,cpp} /tmp/wsclean/msprovider
cp -v parser/*.{h,cpp} /tmp/wsclean/parser
cd /tmp
tar -cjvf ${WORKDIR}/wsclean.tar.bz2 wsclean/
rm -rf /tmp/wsclean
tar -xjvf ${WORKDIR}/wsclean.tar.bz2
cd /tmp/wsclean
mkdir build
cd build
cmake ../
make -j 12
