CC=gcc
CXXC=g++
LIBS=-lm
INCLUDE=-I.
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE

O = cdflib.o

cdflib.a: $(O) 
	ar rcus libcdf.a $(O)

clean:
	rm -f *.o
