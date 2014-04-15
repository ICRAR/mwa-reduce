g++ -o testimpedance -std=c++11 -lcblas -lgsl testimpedance.cpp lnaimpedance.cpp tileimpedance.cpp && ./testimpedance
g++ -o testtile -I/usr/local/include/casacore/ -std=c++11 -lcblas -lgsl -lcfitsio -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lcasa_fits testtile.cpp lnaimpedance.cpp tileimpedance.cpp tilebeam2014.cpp ../fitswriter.cpp -ldl && ./testtile

