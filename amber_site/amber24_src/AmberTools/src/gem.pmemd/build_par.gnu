make clean
make depends
rm config.h
cp config.h.par.gnu config.h
make parallel
