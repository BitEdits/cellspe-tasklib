PROGRAM_ppu   := spunop
IMPORTS        = spu/spunoptest.a -lspe -lc $(SDKLIB)/libsync.a -lpthread
CPPFLAGS = -I/usr/local/include -O4 -g

ifdef CELL_TOP
        include $(CELL_TOP)/make.footer
else
        include $(CELL_SDK_ROOT)/make.footer
endif


