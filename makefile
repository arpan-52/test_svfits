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
CLIBS= -L /usr/lib/x86_64-linux-gnu/ -lz

USE_NOVAS=yes // NOT DEBUGGED DO NOT TURN ON!!!
#DO_DUT1_CORR=yes

ifdef DO_DUT1_CORR
CC += -DUSE_NOVAS -DDO_DUT1_CORR
endif


ifdef USE_NOVAS
CC += -DUSE_NOVAS
STATIC_LIBS= /home/asoft/lib/libcfitsio.a /home/asoft/lib/libnovasc.a
SVOBJ=svfits.o svsubs.o utils.o stats.o novas_prenut.o
else
STATIC_LIBS= /home/asoft/lib/libcfitsio.a /home/asoft/lib/libsla.a
SVOBJ=svfits.o svsubs.o utils.o stats.o sla_prenut.o
endif

ALL= svfits 
all:$(ALL)

svfits:$(SVOBJ) svio.h  gmrt_newcorr.h
	$(CC) $(CFLAGS) -o svfits $(SVOBJ) $(STATIC_LIBS) $(CLIBS)  -lm -lc
clean :
	/bin/rm *.o $(ALL)

