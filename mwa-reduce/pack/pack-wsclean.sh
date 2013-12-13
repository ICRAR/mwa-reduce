WORKDIR=`pwd`
rm -rf /tmp/wsclean
mkdir /tmp/wsclean
mkdir /tmp/wsclean/aocommon
cp CMakeLists-wsclean.txt /tmp/wsclean/CMakeLists.txt
cd ..
cp -v areaset.* banddata.* beamevaluator.* buffered_lane.* cleanalgorithm.* fitsreader.* fitswriter.* imagecoordinates.* imageweights.* inversionalgorithm.* lane.* layeredimager.* multibanddata.* nlplfitter.* matrix2x2.* model.* modelrenderer.* modelsource.* msselection.* polarizationenum.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* stopwatch.* tilebeam.* uvector.* uvwdistribution.* weightmode.* wsclean.* wsinversion.* /tmp/wsclean/
cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h /tmp/wsclean/aocommon
mkdir /tmp/wsclean/parser
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
