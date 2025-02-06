CFLAGS = -g -O3 -Wall -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64

all: vlba14_rearrange

vlba14_rearrange: vlba14_rearrange.o
	$(CC) $(CFLAGS) vlba14_rearrange.c -o vlba14_rearrange

clean:
	rm -f vlba14_rearrange

