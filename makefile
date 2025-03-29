#
# COMPILERS
CC = gcc 
FC = gfortran
#
# FLAGS
#
CFLAGS = -g
#
# Libraries
CLIBS= -L /usr/lib/x86_64-linux-gnu/ -lz
STATIC_LIBS= /home/asoft/lib/libsla.a /home/asoft/lib/libcfitsio.a
#
ALL= svfits 
all:$(ALL)

SVOBJ=svfits.o svsubs.o utils.o stats.o
svfits:$(SVOBJ) svio.h  newcorr.h
	$(CC) $(CFLAGS) -o svfits $(SVOBJ) $(STATIC_LIBS) $(CLIBS)  -lm -lc
clean :
	/bin/rm *.o $(ALL)

