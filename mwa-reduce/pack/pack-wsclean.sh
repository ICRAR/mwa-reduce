WORKDIR=`pwd`
rm -rf /tmp/wsclean
mkdir /tmp/wsclean
cp CMakeLists-wsclean.txt /tmp/wsclean/CMakeLists.txt
cd ..
cp -v areaset.* banddata.* beamevaluator.* cleanalgorithm.* fitswriter.* imagecoordinates.* inversionalgorithm.* lane.* layeredimager.* nlplfitter.* matrix2x2.* model.* modelrenderer.* modelsource.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* tilebeam.* wsclean.* wsinversion.* /tmp/wsclean/
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
