#
# COMPILERS
CC = gcc 
FC = gfortran
#
# FLAGS
#
CFLAGS = -fopenmp -g 
#
# Libraries
#CLIBS= -L /usr/lib/x86_64-linux-gnu/ -lz

#   USE_NOVAS=yes // NOT DEBUGGED DO NOT TURN ON!!!

ifdef USE_NOVAS
CC=gcc -DUSE_NOVAS
STATIC_LIBS= libcfitsio.a libnovasc.a
SVOBJ=svfits.o svsubs.o utils.o stats.o novas_prenut.o
else
STATIC_LIBS= libcfitsio.a libsla.a
SVOBJ=svfits.o svsubs.o utils.o stats.o sla_prenut.o
endif

ALL= svfits 
all:$(ALL)

svfits:$(SVOBJ) svio.h  gmrt_newcorr.h
	$(CC) $(CFLAGS) -o svfits $(SVOBJ) $(STATIC_LIBS) $(CLIBS)  -lm -lc
clean :
	/bin/rm *.o $(ALL)

