# ifndef __NEWCORR_H__
# define __NEWCORR_H__
/* corr.h : Changed on 5 Feb, 1997.
   Fundamental Requirement of structures declared here is that
   they should have same sizes and alignment on all machines.
   DOS is exempted.

   April 1999:  Extensively revised for new version, filename
                changed to newcorr.h
   Mar 2025  :  Multiple other updates which are unfortunately not noted here.
                This file is for use in the SPLOTLIGHT system, derived from
		gvfits-2.05 on 13/Mar/2025. Structures are not changed, but
		some comments are added to make the intent of structure more
		clear, and items that were commented out in the original are
		deleted.
*/

// These are historical labels from the hardware correlator. In later versions
// of the correlator, USB_130 represents the RR polarization and USB_175 the
// LL polarization
enum { USB_130, USB_175, LSB_130, LSB_175, MAX_BANDS };
// Types of stokes products produced RRLL (2 stokes RR and LL) RRRL (4 stokes,
// i.e two copolar and 2 cross polar), RR__ (only one stokes).
enum { RRLL, RRRL, RR__, MACMODES };
// FishCamp Data Acquisition Cards jnc 15/Apr/05 (Historical, unused here,
// except for defining MAX_SAMPS below)
enum { DAS_CARD0, DAS_CARD1, DAS_CARDS, BEDAS_CARDS};
// Historical enum for the Hardware Correlator
enum { IndianPolar, UsbPolar, LsbPolar, UsbCopy, LsbCopy,/* New Mode Names */
       AllU130, AllU175,AllL130, AllL175, 
       arar_arar, alal_alal, brbr_brbr, blbl_blbl,
       aral_brbl, aral_alar, brbl_blbr, arbr_albl,   /* Classic Mode Names */
       DpcMuxVals} ;
enum { NAMELEN=32, DATELEN=32 };
// parameters defining the correlator capacity. There are 32 antennas but only
// 60 samplers (2 per physical GMRT antenna, 2 antennas are dummy). MAC_CARDS
// is historical, right now it is set to MAX_ANTS+1. MAX_BASE appears to be
// 2 times the number of actual distinct baselines (both a1*a2 and a2*a1 are
// counted?)
/* be sure that  MAX_FFTS and MAX_SAMPS are both even numbers ! */
enum { MAX_ANTS=32, MAX_SAMPS=DAS_CARDS*60, MAX_FFTS=MAX_SAMPS, MAC_CARDS=33,  
       MAX_BASE=DAS_CARDS*32*MAC_CARDS, MAX_CHANS=32768, POLS=4  
     };   // %1% MAX_CHANS=256
// parameters relevant for the real time das chain
enum { MAX_PROJECTS=500, MAX_SCANS=512 } ; 
enum { LittleEndian=1, BigEndian=0 };
enum { TransitObs = 32768 } ;
// parameters relevant for the real time das chain, as well as the size of the
// timestamp and weight of each record.
enum { TimeSize=sizeof(double), WtSize = sizeof(double),
       ActiveScans=8, DataFlagSize=ActiveScans*sizeof(int)} ;

// Antenna related parameters including position, fixed delay and phase0.
// samp_id is the sampler-id in the set of samplers currently in use,
// corr->sampler[samp_id].dpc, gives the sampler number in the full master list
// of samplers.
typedef struct
{ char          name[4];
  unsigned char samp_id[MAX_BANDS] ;
  double        bx, by, bz;  /* metres */
  double        d0[MAX_BANDS], p0[MAX_BANDS];
} AntennaParType;

// These are parameters defining the run time correlator setup. These are
// generally defined in the corrsel.def file that is supplied at runtime.
// JNC changed statime to a float (from int for SPOTLIGHT
typedef struct 
{ int    macs, channels, pols, sta, cntrl, iabeam, pabeam1, pabeam2 ;
  float  iabeam_res, pabeam1_res, pabeam2_res, statime,f_step, clock; /* Hz  */
  // clock_used is clock/2<<clk_sel (not used in SPOTLIGHT SVFITS)
  // macmode is the polarization mode, dpcmux,fftmode unused in SVFITS
  unsigned char dpcmux,clksel,fftmode,macmode ;  /* replaces old dummy int */
} CorrParType;


// These are parameters defining the data that is acqtually acquired and
// recorded onto disk. The antmask encodes the antennas for which data
// is acquired (and is available in the SPOTLIGHT raw visibility dump).
typedef struct
{ int   antmask, samplers, baselines, channels, lta, gsb_maxchan, gsb_fstop;
  short bandmask, mode, gsb_stokes;
  short chan_num[MAX_CHANS];  /* %3% Maximum Buffer Allocated */
  /* (mjd_ref*t_unit), (timestamp*t_unit)  in sec */
  double mjd_ref, t_unit, gsb_acq_bw, gsb_final_bw; 
  double dut1;       /* dut1 Assumed from IERS Bulletin D (sec) JNC06Feb09*/
} DasParType ;


// Parameters related to the source being observed. For SPOTLIGHT the antmask
// here would be the one used for forming the real time beams.(?)
typedef struct
{ char object[NAMELEN];
  struct { float i,q,u,v ; } flux ;
  double mjd0 /* fractional mjd, to which ra, dec refer  */ ;
        /*
           mjd0 refers to the epoch of ra_app,dec_app.
           Note that the timestamp in data is wrt to the global
           reference time contained in daspar->mjd_ref
        */
  double ra_app,dec_app,ra_mean,dec_mean,dra,ddec ; /* rad, rad/s */
  double freq[2], first_lo[2],bb_lo[2];   /* Hz */
  double rest_freq[2], lsrvel[2] ;  /* Hz, km/s  */
  double ch_width ;  /* Hz */
  int id, net_sign[MAX_BANDS], mode , dum1;
  unsigned int antmask; 
  unsigned short bandmask, dum2; 
  short calcode, qual ;
} SourceParType ;
// Parameters related to the project that is being run. Antmask here are
// irrelvant for a raw shared memory visibility dump. It is used when forming
// an ltafile using record.
typedef struct
{ char code[8], observer[NAMELEN], title[NAMELEN] ;
  unsigned int antmask ;    /* antennas to record */
  unsigned short bandmask,seq;
} ProjectType ;

// This structure keeps track of the sources and scans that are encountered
// in the course of an observation.
typedef struct
{ int status ;
  float t ;  /* program dependent meaning ! */
  ProjectType proj ;
  SourceParType source ;
} ScanInfoType ;

// Gives the sampler connectivity, i.e. the antenna (ant_id, band) and
// signal chain (fft_id) that a given data stream corresponds to. dpc
// gives the sampler id in the full master list of samplers. The id
// in the SamplerType array in the corr structure below gives the index
// in the array of samplers in use
typedef struct { unsigned char ant_id, band, fft_id, dpc; } SamplerType;

// gives the origin and signal chain of the data streams in each baseline.
// some paramters (card,chip) are historical relating to the hardware
// correlator.
typedef struct
{ SamplerType samp[2] ;
  unsigned char card,chip,pol,word_incr ;
     /*  e.g., RRLL pol=1? word_incr=2 represents RR component ;
               RRLL pol=0?  word_incr=2 represents LL  component ;
               RRRL pol=0? word_incr=2 represents RL component ;
               RR__ pol=0? , word_incr=1
     */
} BaseParType ;         

// This is the master structure containing the full information on the
// correlator, run time setup, data acquisition selections etc.
typedef struct
{ unsigned char endian,dummy[7];
  char          version [NAMELEN   ];  /* should begin with string "CORR" */
  char          bandname[MAX_BANDS ][8];
  AntennaParType  antenna [MAX_ANTS  ];/* Ant pos, freq & other config */
  SamplerType   sampler [MAX_SAMPS ];  /* ant, band vs. ffts           */
  BaseParType   baseline[MAX_BASE  ];  /* Not used in SPOTLIGHT        */
  CorrParType   corrpar;               /* Max. enabled mac_params      */
  DasParType    daspar;                /* Actually useful params       */
} CorrType;

// Rest of this file has various parameters relevant to different modes of
// different correlators that have been use at the GMRT, and are not of
// direct relevance for fits conversion of  spotlight raw visibility dumps
//

// for use in fringe stopping at run time
typedef struct { unsigned int ext_model,idle,stop,userflag ; } AntennaFlagType ;
typedef struct { float phase, dp, delay, dd /* sec/sec */ ; } ModelParType ;
typedef struct { float phase, dp, delay, dslope; } DataParType;
           /* units:  phase in radian, dp in radian/sec,
                      delay in sec, dslope in radians/channel */
typedef struct
{ double t0 ;  /* seconds wrt corr->daspar.mjd_ref */
  int ant_id, band ;
  ModelParType par ;
} ExtModelType ;

# define AntennaTypeSize     sizeof(AntennaParType)
# define MacFftTypeSize      sizeof(MacFftType)
# define SamplerTypeSize     sizeof(SamplerType)
# define BaseParTypeSize     sizeof(BaseParType)
# define DasParTypeSize      sizeof(DasParType)
# define CorrParTypeSize     sizeof(CorrParType)
# define DataParTypeSize     sizeof(DataParType)
# define CorrTypeSize        sizeof(CorrType)
# define CorrSize            sizeof(CorrType)
# define Corr2Size           16192 /* 128 chan, 60 samp CorrType */
# define Corr3Size           29368 /* 256 chan, 120 samp CorrType  ghb*/ 
# define Corr4Size           31152 /* 1024 chan, 120 samp gsb*/ 

#define DEBUG
#undef DEBUG
#define GATHER

#ifdef COMP_BEAM
#define BEAM_MODE1
#define BEAM_MODE2
#endif

#ifdef REDUCED_BANDWIDTH
	#define NTAPS 256 
	#define FIRCHUNKSIZE 1024 
	#define LOWPASS  0
	#define HIGHPASS 1 
	#define BANDPASS 2
	#define BANDSTOP 3
	#define LOWFREQ   0.00
//	#define HIGHFREQ  0.03125
	#define HIGHFREQ  (0.5/DECIMATE)
	#define BARTLETT ippWinBartlett 
	#define BLACKMAN ippWinBlackman
	#define HAMMING  ippWinHamming
	#define HANN     ippWinHann
#endif

#ifndef REDUCED_BANDWIDTH
	#define DECIMATE 1
#endif

#define CHANNEL 256
#define FFTLEN (4*CHANNEL)
#define NCHAN 4
#define CORRLEN (NCHAN*FFTLEN)
#define NUM_POL 16     // Total no of computating nodes // %2% New variable for dual pol
#define NUM_ACQ 16     // No of acquation node          // %2% New variable for dual pol
#define NUM_PROCESSES 32 // Total no of nodes take part
#define NUM_ANT 8        //  No of antennas take part into correlation

#define POLAR_MODE
#ifdef POLAR_MODE
	#define NPOL 2		 // Increase of polar term per baseline	comapre to Intensity mode 
	#define NPOL_T (2*NPOL)  // Total no of polar term per beseline 
	#define NUM_CORR 16      // No of correlation node for a given pol
#else
	#define NPOL 1
	#define NPOL_T 1	
	#define NUM_CORR 8      // No of correlation node for a given pol
#endif

#define FRNG_STEP 0.25
#define NSTEP 2880

#define NCORR (NPOL_T*NUM_ANT*NCHAN*(NUM_ANT*NCHAN+1)/2)
//#define M_PI 3.141592654
#define ACQ_LEN (32*1024*1024)
#define UNPACK_LEN (32*1024*1024)
#define MPI_BUF_CHUNK (ACQ_LEN/NUM_CORR)
#define MPI_EXCESS_CHUNK (64*1024)
#define MPI_OVL_CHUNK (MPI_BUF_CHUNK+MPI_EXCESS_CHUNK)
#define ACQ_OVL (ACQ_LEN+MPI_EXCESS_CHUNK)
#define BEAM_SIZE (UNPACK_LEN/(NCHAN*NUM_ANT))
#define BEAM_SCALE 2*(pow(10,11))

#define SWEEPS_PER_DUMP 1
#define CORR_SIZE ((FFTLEN*NCORR)/2)    //%2% Same size as for one pol!

#define NTCHUNK 16
#define FFTBLOCK 64 

// %%  DAS_BUFSIZE should be corrbuf CORRSIZE.. 
# define DAS_H_KEY 1030
# define DAS_D_KEY 1031
# define DAS_H0_KEY 1032
# define DAS_D0_KEY 1033
// define DAS_BUFSIZE 8192000
# define DAS_BUFSIZE 10240000
# define DAS_HDRSIZE  200000

// %% Required structure..
typedef struct
{ int s0,s1, card;  /* card points to relevant delay/FFT card */
   /* 
      The two relevant samplers are given by 
      daspar.samp_num[s0] and daspar.samp_num[s1]
   */
  int delay, p0, pstep0,p1,pstep1, fstc ;  /* is int ok? */
  /* Do not plant delay difference between two streams here;
     the difference must be handled as static difference
     during initialisation
  */
  float p2fstc, fstc_step ; /* is float fine : sirothia 13oct */
} FftDelayParType ;

typedef struct
{ double clock, t_update;
  double pc_time ;
  double t_off; /* JNC 16Dec98*/
  double delay_step, fstc_scale, nco_tick_time ;
  int cycle, seq, cards ;
            /*cycle = STA cycles between model update */
  unsigned char dpcmux,clksel,fftmode,macmode ;
  ModelParType par[MAX_SAMPS];
  FftDelayParType fdpar[MAX_SAMPS/2];
} ModelInfoType ;

typedef struct
{ int active, status, scan, scan_off;
  CorrType corr ;
  ModelInfoType model ;
  char buf[DAS_HDRSIZE];
} DasHdrType ;

enum { BufMarked=1, BufReady=1<<1,  Rack0Bad=1<<2,Rack1Bad=1<<3,Rack2Bad=1<<4,
       Rack3Bad=1<<5,Rack4Bad=1<<6,Rack5Bad=1<<7, BufFinish=1<<8,
       MaxDataBuf=100
     };

enum { MAX_EVENTS=50000 } ;
typedef struct
{ float t ;
  unsigned char type, cmd ;
  unsigned short seq ;
  int flag_num, scan_num ;  /* indexes on AntennaFlag and ScanInfo */
} EventLogType ;

typedef struct
{ int t_off, wt_off, par_off, data_off, data_words ;
  short par_words, wordsize ;
} RecOffsetType ;

typedef struct
{ int off ;
  BaseParType base ;
  char name[12] ;
} MacWordType ;

typedef struct
{ unsigned char seq, scan ;
  int status,recl, seqnum ;
  RecOffsetType off ;
  MacWordType *mac ;
  float *buf ; /* %1% char *buf */
} RecBufType ;
  
typedef struct
{ int active,status;
  unsigned short events,flags,starts,stops ;
  CorrType corr ;
  AntennaFlagType flag[MAX_EVENTS][MAX_BANDS] ;
  EventLogType event[MAX_EVENTS] ;
  ScanInfoType scaninfo[MAX_SCANS] ;
  RecOffsetType offset ;
} DataInfoType ;

/* GSBE : DatabType, DataBufType not required structures */
typedef struct
{ int flag,rec,seqnum ;
  unsigned short flag_seq, newstate ;
} DataTabType;

typedef struct
{ int flag, blocksize, maxblocks, cur_block, first_block, cur_rec;
  DataTabType dtab[MaxDataBuf];
// char buf[DAS_BUFSIZE] ; //
   char buf[DAS_BUFSIZE] ;      
} DataBufType ;

# endif
