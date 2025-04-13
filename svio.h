#ifndef __SVIO_H__ 
#define __SVIO_H__ 

#include <stdio.h>

#ifndef __NEWCORR_H__ 
#include "newcorr.h"
#endif

// added this for SPOTLIGHT, not clear where it was getting defined earlier
// assume that this code will run only on 64 architecture
#define ARCH_64BIT  
#define PATHLEN 256
#define LINELEN  80

typedef struct { float r, i ; } Complex ;
typedef struct { float r, i, wt ; } Cmplx3Type ;
typedef struct { double bx,by,bz, u,v,w ; } BaseUvwType ;

// Generic Structure defining the details of some given random parameter
typedef struct{
/* 
   standard "type"s in RandParType :
   A=ascii, B=byte=unsigned char, S=short, I=int,
   L=long,  F=float, D=double C=Complex {float real,imag;},

   Extended "type"s under consideration:
   c2=C={float real,imag ;}
   c3={float real,imag,weight;}, an=AntennaType, co=CorrType,
   sc=ScanParType, cp=CorrParType, dt=DataParType, ba=BaseType
   da=DasParType, sa=SamplerType, cn=ConnectionType
*/
  char code[8] ;
  int form,off, nmemb ;
  /* off = byte offset (zero-relative) in the buffer */
  /* "nmemb" parameters, whose element size is determined by form  */
} RandParType ;

// uvw and other parameters that are associated with a FITS visibility record
typedef struct { float u,v,w, baseline, date1,date2, su, fq ; } UvwParType ;

typedef struct { 
  unsigned int off;  // offset to start of data for this baseline
  int          drop; // flagged?
  char         ant0, ant1,band0,band1;
  int          flip; // flip the imaginary part to match AIPS/CASA convention
} VisInfoType ;
// order of baselines in the output FITS file
typedef struct vis_par_type 
{ int         vis_sets; // numBaselines to combine to 1 FITS group (==numSTOKES)
  int         blocks ;  // number of records to process at a time  (?)
  int         if_stride, chan_stride, pol_stride,block_stride ;//unused
  int         pol_ind[8],if_ind[8], set_off[8] ;//unused
  struct      {int pol0,pol1 ; } self_pol[8] ;//unused
  VisInfoType visinfo[MAX_BASE] ;
} VisParType ;
typedef struct{
  int dataform, samplers,channels, nself,ncross ;
  RecOffsetType off ;
  MacWordType word[MAX_BASE] ;
} MacInfoType ;

enum{RawPars=1}; //ADHOC CHECK IF RAWPARS IS USED - SEEMS FLAG RELATED
typedef struct init_hdr_type{
  int  file, scans, type, valid, recl, hdr_recs ;
  RecOffsetType off ;
  RandParType item[RawPars] ; /* valid items defined among these */
  CorrType corr ;
  char rec_form[NAMELEN], obs_mode[NAMELEN], version[NAMELEN];
  MacInfoType mac;
} InitHdrType;

typedef struct scan_rec_type {
  double t_start,t_end ;
  unsigned char array_id,scan_id, source_id,freq_id ;
  int file, scannum, hdr_rec, records, firstrec ,lastrec ;
  CorrType *corr ;
  ScanInfoType *scan ;    /* may2000  */
} ScanRecType ;

enum{MaxRecFiles=16};
typedef struct rec_file_par_type{
  /* information of the raw timesliced data files*/
  int     nfiles;              // total number of raw data files
  char    path[PATHLEN];       // path to the raw data files
  char    fname[MaxRecFiles][LINELEN]; // names of the raw data files
  double  t_start[MaxRecFiles];// start time (IST?) of each file
  double  b_start[MaxRecFiles];// start time (IST?) of first rec with data
  double  t_slice;             // slice duration (sec)
  double  slice_interval;      // interval between 2 slices in file (sec)
  int     rec_per_slice;       // number of lta records per slice
} RecFileParType;
  
typedef struct burst_par_type{
  char   name[LINELEN];
  double mjd;      // burst time (mjd)
  char   date[LINELEN];//burst date
  double t,dt;     // burst time, error in burst time (sec) [IST]
  double width;    // observed burst width (sec) [Dispersion dominated]
  double int_wd;   // intrinsic width of the burst (sec)
  double DM,dDM;   // DM and error in DM
  double f;        // ref. freq for burst time (Hz)
  int    bm_id;    // id of the beam containing the burst
  double ra_app,dec_app;   // apparent ra,dec of the beam containing the burst;
} BurstParType;

typedef struct bpass_type{
  Cmplx3Type *mean[MAX_BASE];// mean of re,im over the slice
  float      *abp[MAX_BASE]; // normalized amplitude of the bandpass
}BpassType;
typedef struct sv_selection_type
{ int              scans, baselines , antennas ;
  unsigned int     antmask ;
  short            stokes,sidebands,channels,sideband[2], stokes_allow[4] ;
  short            start_chan,force_app,chan_inc;
  short            stokes_type,fake_data,dum[2];
  char             fitsfile[LINELEN] ;
  float            iatutc,fdum ;
  double           epoch,timestamp_off;
  char             coord_type[NAMELEN];
  VisParType       vispar ;
  BpassType        bpass;
  InitHdrType     *hdr  ;
  CorrType        *corr ;
  ScanRecType     *srec ;
  int              do_log;
  FILE            *lfp;
  RecFileParType   recfile;  
  BurstParType     burst;
  int              update_burst; //compute burst params (see update_burst())
  int              do_flag; // mad based flagging of data 
  float            thresh; // threshold for flagging (units of mad)
  int              all_chan;//ignore BurstPar and copy all chans (debug use)
  int              do_band;// apply amplitude bandpass calibration

} SvSelectionType ;

// various structures used while making FITS file tables
#pragma pack(push,4) /* push current alignment to stack */
struct an_struct
       { char  anname[8] ;
         double x,y,z ;
#ifdef ARCH_64BIT
         int nosta,mntsta ;
#else
         long nosta,mntsta ;
#endif
         float staxof ;
         char poltya[4]  ;
         float polaa, polcala[4] ;
         char poltyb[4] ;
         float polab, polcalb[4] ;
       } ;
#pragma pack() /* push current alignment to stack */

#pragma pack(push,4) /* push current alignment to stack */
struct su_struct
       { int id ;
         char source[16] ;
         int qual ;
         char calcode[4] ;
         float iflux[2],qflux[2],uflux[2],vflux[2];
         double freq_off[2],bandwidth, ra,dec,epoch,ra_app,dec_app ;
//         double lsrvel[2], rest_freq[2], pmra,pmdec ;
         double lsrvel, rest_freq, pmra,pmdec ; //JNC8apr16 for Casa4.2.5
       } ;
#pragma pack() /* push current alignment to stack */
                  
#pragma pack(push,4) /* push current alignment to stack */
struct fq_struct 
       { int id ;
         double   iffreq[1];
         float    chwidth[1];
         float    bandwidth[1];
#ifdef ARCH_64BIT
         int sideband[1];
#else
         long int sideband[1];
#endif
       } ;
#pragma pack() /* push current alignment to stack */

#pragma pack(push,4) /* push current alignment to stack */
struct fq_struct2 
       { int id ;
         double   iffreq[2];
         float    chwidth[2];
         float    bandwidth[2];
//         long int sideband[2];
         int sideband[2];
       } ;
#pragma pack()  /* restore original alignment from stack */

typedef struct
{ int npar,maxpoints ;
  double t_int, t_ref ;
  Complex *databuf;
  double *tbuf ;
  float *rbuf,*ibuf, *xbuf, *parbuf ;
} FitBufType ;

int copy_all_chans(SvSelectionType *user, int idx, char **outbuf, int restart);
int     copy_sel_chans(SvSelectionType *user, int idx, char **buf,int *restart);
int     clip(char *visbuf, SvSelectionType *user, int groups);
int     fake_sel_chans(SvSelectionType *user, int idx, char **buf,int *restart);
unsigned short  float_to_half(const float x);
float   half_to_float(const unsigned short x);
int     init_user(SvSelectionType *user, char *fname, char *fname1);
double  lmst(double mjd);
char   *mjd2iau_date(double mjd);
float   svVersion(void);
void    swap_bytes(unsigned short *, int) ;
void    swap_short(unsigned *, int) ;
void    swap_long(void *, int) ;
void    swap_d(double *, int) ;
int     robust_stats(int n, float *x, float *med, float *mad);


#endif
