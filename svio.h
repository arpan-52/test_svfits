#ifndef __SVIO_H__ 
#define __SVIO_H__ 

#include <stdio.h>

#ifndef __NEWCORR_H__ 
#include "gmrt_newcorr.h"
#endif

// added this for SPOTLIGHT, not clear where it was getting defined earlier
// assume that this code will run only on 64 architecture
#define ARCH_64BIT  
#define PATHLEN  1024
#define LINELEN  1024
#define MAX_REC_PER_SLICE 50  // records in a given slice

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
  int     have_idx;//          // format includes index if true
  char    path[PATHLEN];       // path to the raw data files
  char    fname[MaxRecFiles][LINELEN]; // names of the raw data files
  double  t_start[MaxRecFiles];// start time (IST?) of each file
  double  b_start[MaxRecFiles];// start time (IST?) of first rec with data
  int     start_rec[MaxRecFiles];// first record with burst signal
  int     n_rec[MaxRecFiles];  // number of contiguous records with burst signal
  double  t_slice;             // slice duration (sec)
  double  slice_interval;      // interval between 2 slices in file (sec)
  int     rec_per_slice;       // number of lta records per slice
  double  mjd_ref;             // the reference mjd for the raw visibility file
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
  int         file_idx,slice;   // file index and slice of this data
  float      *abp[MAX_BASE];    // normalized amplitude of the "bandpass"
  Complex    *off_src[MAX_BASE];// mean of re,im over off source region
  short       start_chan[MAX_REC_PER_SLICE];// per rec start channel of burst
  short       end_chan[MAX_REC_PER_SLICE];  // per rec end channel of signal 
}BpassType;

typedef struct dut1_table_type { double mjd; double dut1;} DUT1TabType;
enum{DUT1_TABSIZE=1024};
#define DUT1_PREDICTION_BUFFER  7.0 /* WARNING ISSUED IF PREDICTION AVAILABLE*/
                                    /* FOR FEWER THAN THESE MANY DAYS       */
typedef struct sv_selection_type
{ int              scans, baselines , antennas ;
  unsigned int     antmask ;
  short            stokes,sidebands,channels,sideband[2], stokes_allow[4] ;
  short            start_chan,force_app,chan_inc;
  short            stokes_type,fake_data,dum[2];
  float            statime; //statime in corr is int and not suitable
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
  FILE            *lfp,*pcfp;
  RecFileParType   recfile;  
  BurstParType     burst;
  int              update_burst; //compute burst params (see update_burst())
  int              do_flag; // mad based flagging of data 
  float            thresh; // threshold for flagging (units of mad)
  int              all_chan;//copy all chans for records with burst
  int              all_data;//copy all data; ignore burstpar
  int              nchav; // number of channels to average (all_data only)
  int              do_band;// apply amplitude bandpass calibration
  int              do_base;// remove baseline (i.e. mean offsource visibilities)
  int              postcorr;//generate postcorr beam output
  int              drop_csq;// drop CSQ baselines while making postcorr beam
  int              num_threads;//number of threads to use in parallel sections
  double           i2eapp[3][3],emean2i[3][3];//matrices for J2000 coversion
  int              recentre; // recentre visibilities to new phase centre
  double           rmat[3][3]; //rotation matrix for new phase centre
  double           lmn[3]; // offset coordinates - original UVW frame(J2000)
  double           lmn_a[3]; //offset coordinates - original uvw frame (app)
  DUT1TabType      dtab[DUT1_TABSIZE];
  int              n_dut1;
  double  lta;                 // lta time (sec) (only when converting all data)
  int     n_lta  ;             // num of output lta records (-1==> all)
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

typedef unsigned short ushort;

int     copy_vis(SvSelectionType *user, BpassType *bpass, int idx, int slice,
		 int start_rec, int n_rec,char *rbuf,char **outbuf);
int     avg_vis(SvSelectionType *user, int idx, int slice, char *rbuf,
		char *outbuf);
int     clip(char *visbuf, SvSelectionType *user, int idx, int slice,
	     int groups);
int     fake_sel_chans(SvSelectionType *user, int idx, char **buf,int *restart);
ushort  float_to_half(const float x);
int     get_file_order(SvSelectionType *user, int *order);
double  get_ha(SvSelectionType *user, double tm);
float   half_to_float(const unsigned short x);
int     init_dut1tab(SvSelectionType *user, char *bulletinA);
void    init_mat(SvSelectionType *user, double tm);
int     init_user(SvSelectionType *user, char *fname, char *anthdr,
		  char *bhdrfile, char *bulletinA);
double  lmst(double mjd);
char   *mjd2iau_date(double mjd);
int     read_slice(SvSelectionType *user, int idx, int slice, char *rbuf);
float   svVersion(void);
void    swap_bytes(unsigned short *, int) ;
void    swap_short(unsigned *, int) ;
void    swap_long(void *, int) ;
void    swap_d(double *, int) ;
int     robust_stats(int n, float *x, float *med, float *mad);

#ifdef USE_NOVAS
void novas_app2j2000(double rap, double decp, double mjd, double iatutc,
	       double *ra,double *dec);
int  novas_mean2j2000(double rap, double decp, double mjd, double iatutc,
	       double *ra,double *dec);
int  novas_prenut_vis(SvSelectionType *user,UvwParType *uvw, double mjd);
#else
void   sla_amp_(double *ra_app, double *dec_app, double *mjd, 
	      double *epoch1,double *ra_mean, double *dec_mean);
void   sla_preces_(char *system,double *epoch, double *epoch1,
		 double *ra, double *dec,int len);
void sla_prenut_vis(UvwParType *uvw,double mjd,double ra_app,double dec_app,
		    double epoch1);
void   sla_nut_(double *mjd, double *a); 
double sla_epj_(double *mjd);

#endif //USE_NOVAS

#endif
