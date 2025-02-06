
# PPU_COMPILER = xlc

# PROGRAM_ppu64   := minitest-cell
PROGRAM_ppu     := minitest-cell

IMPORTS         = cellspe-tasklib.a -lspe -lc $(SDKLIB)/libsync.a -lfftw3f -lpthread
CPPFLAGS = -I/usr/local/include -O4 -g -pthread -fpermissive -DUSE_CELL_SPE=1 -m32 -maltivec

# for debug code: (-O0), try to avoid -pg profiling info generation on PS3
# CPPFLAGS = -g -pg -I/usr/local/include/
# LDFLAGS  = -pg

ifdef CELL_TOP
        include $(CELL_TOP)/make.footer
else
        include $(CELL_SDK_ROOT)/make.footer
endif


