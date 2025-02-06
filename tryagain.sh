
# recompiles SPE code, PPU code, executes test sets 

make -f Makefile-cellspetasklib clean
make -f Makefile-cellspetasklib ; rm cellspe-tasklib.o  ; make
./minitest-cell
