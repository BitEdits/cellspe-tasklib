#!/bin/sh

# This script compiles the matrix multiplication example.

# change this according to your needs
# SDK 2.1, 2.0 or directly on Cell system
# CELL_BIN="/opt/cell/bin"
# CELL_BIN="/opt/cell/toolchain-3.3/bin"
CELL_BIN="/usr/bin"
SPE_INCLUDES="-I/opt/ibm/cell-sdk/prototype/sysroot/usr/spu/include/"
SPE_LIBS="-L/opt/ibm/cell-sdk/prototype/sysroot/usr/spu/lib/ -lgmath -lsimdmath -lmisc"

# remove previously compiled binary
rm -f minicorrelator minicorrelator_spu

# compile SPE code
echo "${CELL_BIN}/spu-gcc -W -Wall $SPE_INCLUDES -include spu_intrinsics.h -O3 -g -c minicorrelator_spu.c ${SPE_LIBS}"
${CELL_BIN}/spu-gcc -W -Wall $SPE_INCLUDES -include spu_intrinsics.h -O3 -g -c minicorrelator_spu.c ${SPE_LIBS}

 # -S -dA
 # /opt/ibm/cell-sdk/prototype/bin/spu_timing minicorrelator_spu.s

# link SPE .o files together
# echo "${CELL_BIN}/spu-gcc -o minicorrelator_spu minicorrelator_spu.o minicorrelator_spu_simd.o"
# ${CELL_BIN}/spu-gcc -o minicorrelator_spu minicorrelator_spu.o minicorrelator_spu_simd.o
echo "${CELL_BIN}/spu-gcc -o minicorrelator_spu minicorrelator_spu.o ${SPE_LIBS}"
${CELL_BIN}/spu-gcc -o minicorrelator_spu minicorrelator_spu.o ${SPE_LIBS}

# compile and link SPE eval code
# ${CELL_BIN}/spu-gcc -W -Wall $SPE_INCLUDES -include spu_intrinsics.h -O3 -g -c eval_spu.c ${SPE_LIBS}
# ${CELL_BIN}/spu-gcc -o eval_spu eval_spu.o ${SPE_LIBS}

# embedd SPE object file into PPE object
echo "${CELL_BIN}/ppu-embedspu -m64 minicorrelator_spu minicorrelator_spu minicorrelator_spu-embed64.o"
${CELL_BIN}/ppu-embedspu -m64 minicorrelator_spu minicorrelator_spu minicorrelator_spu-embed64.o

# combile PPE code
echo "${CELL_BIN}/ppu-gcc -W -Wall -O3 -g -I /opt/ibm/cell-sdk/prototype/src/lib/misc/ -c minicorrelator_ppu.c"
${CELL_BIN}/ppu-gcc -W -Wall -O3 -g -I /opt/ibm/cell-sdk/prototype/src/lib/misc/ -c minicorrelator_ppu.c

# link SPE adn PPE object files together
echo "${CELL_BIN}/ppu-gcc -o minicorrelator minicorrelator_ppu.o minicorrelator_spu-embed64.o -lspe"
${CELL_BIN}/ppu-gcc -o minicorrelator minicorrelator_ppu.o minicorrelator_spu-embed64.o -lspe

# rm -f minicorrelator_spu *.o
rm -f *.o
