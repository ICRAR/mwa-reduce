g++ -o testimpedance -std=c++11 -lcblas -lgsl -lcfitsio -lboost_filesystem -lboost_system testimpedance.cpp lnaimpedance.cpp tileimpedance.cpp system.cpp ../fitsiochecker.cpp && ./testimpedance
g++ -o testtile -I/usr/local/include/casacore/ -std=c++11 -lcblas -lgsl -lcfitsio -lcasa_ms -lcasa_tables -lcasa_casa -lcasa_measures -lcasa_fits -lboost_filesystem -lboost_system testtile.cpp lnaimpedance.cpp tileimpedance.cpp tilebeam2014.cpp joneslookupdipole.cpp system.cpp ../fitswriter.cpp ../fitsiochecker.cpp ../alglib/*.cpp -ldl && ./testtile

