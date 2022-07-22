CC=gcc
CXXC=g++
LIBS=-lm
INCLUDE=-I.
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE

O =  energy_par.o fold.o fold_vars.o utils.o params.o gquad.o aln_util.o duplex.o cofold.o subopt.o alifold.o list.o ribo.o

libRNAfold.a: $(O) 
	ar rcus libRNAfold.a $(O)
	ranlib libRNAfold.a

energy_par.o: energy_par.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c energy_par.c

fold.o: fold.c fold.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c fold.c

fold_vars.o: fold_vars.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c fold_vars.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c utils.c

gquad.o: gquad.c gquad.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c gquad.c

aln_util.o: aln_util.c aln_util.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c aln_util.c
	
params.o: params.c
	$(CC) $(CFLAGS) $(HG_DEFS) -c params.c

duplex.o: duplex.c duplex.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c duplex.c
	
cofold.o: cofold.c cofold.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c cofold.c

subopt.o: subopt.c subopt.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c subopt.c
	
alifold.o: alifold.c alifold.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c alifold.c

list.o: list.c list.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c list.c

ribo.o: ribo.c list.h
	$(CC) $(CFLAGS) $(HG_DEFS) -c ribo.c
	
clean:
	rm -f *.o libRNAfold.a


