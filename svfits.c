/*
  Main program for converting SPOTLIGHT raw visibility dump files to a random
  group UV FITS file. Derived from gvfits.c

  The first few functions creates the header and the tables associated with
  the FITS file (tables contain the details of the antenna locations,
  source parameters, frequency settings etc.). Some default values for these
  parameters are hardcoded in the program, but the intention is that they
  would be read from meta data files created by the SPOTLIGHT pipeline at
  the time of dumping the visibilities.

  The main program first sets up the structures holding the meta data (these
  are essentially the structures in newcorr.h, which are also used by the
  real time programs), as well as the structures holding information required
  for the format conversion, and the parameters of the burst for which the
  visibility data has been dumped, then coverts the visibilities into random
  groups, and finally creates the associated tables.

  Only the visibilities containing the burst signal are extracted, for more
  details see copy_sel_chans().

  JNC March 2025

  added code to limit the amount of memory used at a given time, files can
  be processed chunk by chunk in case the memory needed to process the entire
  burst data in a given file exceeds the limit. (This may not be really needed,
  and I might remove it in future to keep the code from getting unneccessarily
  complex).

  added the option to copy all channels, instead of just the channels containing
  the burst signal. This produces a regular multi-channel random groups UV file.
  This mode is helpful mainly for debuging. In this mode options to do
  amplitude bandpass calibration, and to remove the mean visibility over
  the entire slice have also been provided. The idea is to provide them also
  in the regular mode (i.e. where only the burst signal is copied into a
  single channel output file). Also added the option of updating the burst
  parameters based on the burst MJD and the burst intrinsic width (these
  are the two parameters expected from the SPOTLIGHT pipeline).

  jnc apr 2025

  major clean up of the code, the major difference is that it now processes
  a "slice" (i.e. the 50 continuguous records in any given raw visibilty file)
  at a time. This makes is more straightforward to do amplitude bandpass
  calibration, subtract RFI etc. It also helps make the code much more modular,
  allowing for multiple operations to be done on the slice once the slice has
  been read into memory.
  
  An option to produce a post-correlation beam (currently at the phase centre,
  in a later version the plan is to rotate it to the centre of the realtime
  beam in which the burst was located) has also been added, along with various
  other options. It is easiest to look at the associated sv_par.fits file to
  understand the options. The log file has also been reformatted to make it
  easier for e.g. to overplot the data actually selected for the burst on top
  of the post_correlation beam data produced by this program. Currently it is
  a little slow, it would probably be good to parallelize some of the visibility
  handling.

  jnc apr 2025
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <endian.h>
#include <unistd.h>

#include "fitsio.h"
#include "fitsio.h"
#include "longnam.h"
#include "newcorr.h"
#include "svio.h"

static float RefFreq,ChanWidth;


int replace_nulls(char *str, int len)
{ int i,k ;
  k = strlen(str) ;
  if (k > len) return 0 ;
  for (i=k; i<len; i++) str[i] = ' ' ;
  return i ;
}

double jd2gst(double jd)  // gst in degrees corresponding to jd
{
  double c0,c1,c2,c3;
  double temp,t;
  c0=24110.54841;
  c1=8640184.812866;
  c2=0.093104;
  c3=-6.2E-6;
  t=(jd-2451545.0)/36525;
  temp=((c0+c1*t+c2*t*t+c3*t*t*t))/240.0;
  return (temp);
}
/* create the antenna table
   jnc mar25 - made some minor changes to avoid the overflow warning from gcc
*/
int writeagtable(fitsfile *fptr, SvSelectionType *user, float iatutc)
{ double gstia0 ;
  ScanRecType *srec = user->srec ;
  CorrType *corr = user->corr ;
  char      tstr[256];
  int row, status;
  int tfields =12;
  long nrows;
  char ref_date[24] ;
  char *ttype[]={ "ANNAME","STABXYZ","ORBPARM","NOSTA","MNTSTA","STAXOF",
		    "POLTYA","POLAA","POLCALA","POLTYB","POLAB","POLCALB"};
  char *tform[]={"8a","3d","0d","1j","1j","1e","4a","1e","4e","4a","1e","4e"};
  char *tunit[]={" ", "METERS", " "," "," ","METERS", " ", "DEGREES", " ", " ",
		   "DEGREES", " "};

  struct an_struct agt ;

  status=0;
  strncpy(ref_date , mjd2iau_date(corr->daspar.mjd_ref), 10) ;
  ref_date[10] = 0 ;

  { double mjd_iat0 = corr->daspar.mjd_ref +5.5/24.0 + iatutc/3600.0 ;
    gstia0 = jd2gst(mjd_iat0+2400000.5) ;
  }

  if(strstr(corr->version,"COR30x2 ")==NULL && strstr(corr->version,"HST03") ==NULL && strstr(corr->version,"HST04") ==NULL ) {
    nrows=30; // MAX_ANTS of old corr Pre_POLAR
  }
  if(strstr(corr->version,"COR30x2 ")!=NULL && strstr(corr->version,"HST02") !=NULL) {
    nrows=30; // MAX_ANTS of old corr POLAR
  }
  if(strstr(corr->version,"COR30x2 ")!=NULL && strstr(corr->version,"HST03") ==NULL) {
    nrows=30; // MAX_ANTS of old corr POLAR
  }
  if(strstr(corr->version,"COR30x2 ")!=NULL && strstr(corr->version,"HST03") !=NULL) {
    nrows=MAX_ANTS; // MAX_ANTS of GSB = 32
  }
  if(strstr(corr->version,"GWB-III ")!=NULL && strstr(corr->version,"HST04") !=NULL) {
    nrows=MAX_ANTS; // MAX_ANTS of GWB = 32
  }
  if(strstr(corr->version,"SPOTLIGHT")!=NULL){
    nrows=MAX_ANTS; // MAX_ANTS of SPOTLIGHT = 32
  }


  if ( ffcrhd(fptr,&status)) return status ;
  if (ffphbn (fptr, nrows,tfields,ttype,tform,tunit,"AIPS AN",0L,&status))
    return status ;
  if (ffpkyj(fptr,"EXTVER", 1,"Version number of table",&status))
    return status ;
  if (ffpkyd(fptr,"ARRAYX", 1657004.6290,12,"",&status))
    return status ;
  if (ffpkyd(fptr,"ARRAYY", 5797894.3801,12,"",&status))
    return status ;
  if (ffpkyd(fptr,"ARRAYZ", 2073303.1705,12,"",&status))
    return status ;
  if (ffpkyd(fptr,"GSTIA0", gstia0,12,"GST at IAT=0 (degrees) on ref. date",&status))
    return status ;
  if (ffpkyd(fptr,"DEGPDY",0.36098564497436661e3 ,20,
	     "Earth rotation rate (deg/IAT day",&status))
    return status ;
  if (ffpkyd(fptr,"FREQ",RefFreq,12,"Reference frequency for array",&status))
    return status ;
  if (ffpkys(fptr,"RDATE", ref_date,"Reference date 'YYYY-MM-DD'",&status))
    return status ;
  if (ffpkyd(fptr,"POLARX", 0.0,12,"Polar position X (meters) on Ref. date",&status))
    return status ;
  if (ffpkyd(fptr,"POLARY", 0.0,12,"Polar position Y (meters) on Ref. date",&status))
    return status ;
  if (ffpkyd(fptr,"UT1UTC", 0.0,12,"UT1-UTC (time sec.)",&status))
    return status ;
  if (ffpkyd(fptr,"DATUTC", 0.0,12," ",&status))
    return status ;
  if (ffpkys(fptr,"TIMSYS", "IAT","Time system, 'IAT' or 'UTC' ",&status))
    return status ;
  if (ffpkys(fptr,"ARRNAM", "GMRT","Array name",&status))
    return status ;
  if (ffpkyj(fptr,"NUMORB", 0,"Number of orbital parameters",&status))
    return status ;
  if (ffpkyj(fptr,"NOPCAL", 4,"Number of pol. cal constants",&status))
    return status ;
  if (ffpkyj(fptr,"FREQID",-1 ,"The ref freq of the uv data",&status))
    return status ;
  if (ffpkyd(fptr,"IATUTC", iatutc,12,"IAT-UTC (time sec) ",&status))
    return status ;
  if (ffpkys(fptr,"POLTYPE", "APPROX","Feed polarization parameterization" ,&status))
    return status ;
  if (ffpkyj(fptr,"P_REFANT" ,5, "Reference antenna" ,&status))
    return status ;
  if (ffpkyd(fptr,"P_DIFF01", 0.0,12," ???",&status))
    return status ;
  if (ffpkyd(fptr,"P_DIFF02", 0.0,12," ???",&status))
    return status ;

  /* initialize and write the binary table */
  for(row=1; row <= nrows; row++ )
  { AntennaParType *ant = user->corr->antenna + row - 1;
    strncpy(agt.anname,ant->name,3);agt.anname[3]='\0';
    sprintf(tstr,":%02d", row) ;
    strncat(agt.anname,tstr,3);
    agt.anname[6] = agt.anname[7] = ' ' ;
    agt.x = ant->bx;
    agt.y = ant->by;
    agt.z = ant->bz; 
    agt.nosta=row;
    agt.mntsta=0;
    agt.staxof=0.0;

    if (srec->scan->source.freq[0] < 850e6)
            { agt.poltya[0] = 'R' ;   agt.poltyb[0] = 'L'; }
    else    { agt.poltya[0] = 'X' ;   agt.poltyb[0] = 'Y'; }
    agt.poltya[1] = agt.poltyb[1] = ' ' ;
    agt.poltya[2] = agt.poltyb[2] = ' ' ;
    agt.poltya[3] = agt.poltyb[3] = ' ' ;

    agt.polaa=0.0;
    agt.polcala[0]=0.0;
    agt.polcala[1]=0.0;
    agt.polcala[2]=0.0;
    agt.polcala[3]=0.0;

    agt.polab=0.0;
    agt.polcalb[0]=0.0;
    agt.polcalb[1]=0.0;
    agt.polcalb[2]=0.0;
    agt.polcalb[3]=0.0;
    if(corr->endian & LittleEndian)
    {
      swap_d(&agt.x,3) ;
      swap_long(&agt.nosta, 3) ;
      swap_long(&agt.polaa, 5) ;
      swap_long(&agt.polab, 5) ;
    }
//	fprintf(stdout, "SIZE OF AN = %d\n", sizeof(struct an_struct));
   	if (ffptbb(fptr, row, 1, sizeof(struct an_struct), (unsigned char *) &agt, &status)) return status ;
  }

  return status ;
}
// create the frequency table
int writefqtable(fitsfile *fptr, SvSelectionType *user, int scans, int frequencies)
{
  enum {tfields = 5 } ;
  CorrType *corr = user->corr ;
  ScanRecType *srec = user->srec ;
  SourceParType *source ;
  int status ;
  int i,j, row ;
  struct fq_struct fqt ;
  struct fq_struct2 fqt2 ;

  char *ttype[tfields] =
     { "FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND" };
  char *tform[tfields] = 
     {   "1J" ,    "1D" ,      "1E" ,          "1E" ,          "1J" } ;
  char *tform2[tfields] = 
     {   "1J" ,    "2D" ,      "2E" ,          "2E" ,          "2J" } ;
  char *tunit[tfields] = {"  ", "Hz", "Hz", "Hz" , "  " } ;

  status=0;

//  fprintf(stderr, "### No. of Freq = %d\n", frequencies);
  if ( ffcrhd(fptr,&status)) return status ;
  if(user->sidebands==1)
  { if (ffphbn (fptr, frequencies,tfields,ttype,tform,tunit,
		"AIPS FQ", 0L, &status)) return status ;}
  else
  { if(user->sidebands==2)
    { if (ffphbn (fptr, frequencies,tfields,ttype,tform2,tunit,
		  "AIPS FQ", 0L, &status)) return status ;
    }
    else{ fprintf(stderr,"Illegal number of sideband %d\n",user->sidebands); return -1;}
  }
  if (ffpkyj(fptr,"EXTVER", 1,"Version number of table",&status))
    return status ;
  if (ffpkyj(fptr,"NO_IF" ,user->sidebands, "The number of IFs" ,&status))
    return status ;

  for (row=1; row <= frequencies; row++)
  {
    for (i=0; i<scans; i++)if(srec[i].freq_id == row) break ;
    if (i == scans) 
      { fprintf(stderr,"FATAL: could not locate freq %d\n", row) ;
      return -1 ;
    }

    source = &srec[i].scan->source ;
    fqt.id = fqt2.id=srec[i].freq_id ;
//  fprintf(stdout, "FREQ ID  ..... = %d\n", fqt.id);

    for (j=0; j<user->sidebands; j++) 
    { int sideband=corr->baseline[0].samp[0].band/2;  /* 0=>USB, 1=>LSB */
      int pol=user->stokes_allow[0]/2;  /* 0=>130, 1=>175 */
      int sc=corr->daspar.chan_num[user->start_chan],band;
//      fprintf(stdout, "SIDEBANDS..... = %d\n", user->sidebands);
      if(user->sidebands==1)
      { band=2*sideband+pol;
	fqt.iffreq[j] = source->freq[pol] + sc*source->ch_width*source->net_sign[band] - RefFreq ;
        fqt.chwidth[j] = user->chan_inc*source->ch_width*source->net_sign[band]; /* negative for lsb*/
	fqt.bandwidth[j] = (user->channels*user->chan_inc*source->ch_width);
	fqt.sideband[j] = 1 ; /*always usb */

      }else
      {  band=2*user->sideband[j]+pol;
	fqt2.iffreq[j] = source->freq[pol] + sc*source->ch_width*source->net_sign[band] - RefFreq ;
        fqt2.chwidth[j] = user->chan_inc*source->ch_width*source->net_sign[band];/*negative for lsb*/
	fqt2.bandwidth[j] = user->channels*user->chan_inc*source->ch_width ;
	fqt2.sideband[j] = 1 ;  /* always usb*/
      }
    }
    if (corr->endian & LittleEndian)
    { 
      if(user->sidebands==1)
      { swap_long(&fqt.id, 1) ;
        swap_d(fqt.iffreq, user->sidebands) ;
	swap_long(fqt.chwidth,3*user->sidebands) ;
      }else
      { swap_long(&fqt2.id, 1) ;
        swap_d(fqt2.iffreq, user->sidebands) ;
	swap_long(fqt2.chwidth,3*user->sidebands) ;
      }
    }

    if(user->sidebands==1) {
      ffptbb(fptr, row,1, sizeof(struct fq_struct), (unsigned char *)&fqt, &status) ;
    }
    else
      ffptbb(fptr, row,1, sizeof(struct fq_struct2), (unsigned char *)&fqt2, &status) ;
  }
  return status ;
}
// create the source table
int writesutable(fitsfile *fptr, SvSelectionType *user, int scans, int sources, double epoch1) 
{
  enum {tfields=19};
  ScanRecType *srec = user->srec ;
  SourceParType *source ;
  BurstParType  *burst=&user->burst;
  CorrType *corr = user->corr ;
  struct su_struct sut ;
  double ra,dec,epoch ;
  int status ;
  int i,j, row ;
  char *ttype[tfields] = 
       { "ID. NO. ", "SOURCE  ", "QUAL    ", "CALCODE ", "IFLUX   ",
         "QFLUX   ", "UFLUX   ", "VFLUX   ", "FREQOFF ", "BANDWIDTH",
         "RAEPO   ", "DECEPO  ", "EPOCH   ", "RAAPP   ", "DECAPP  ",
         "LSRVEL  ", "RESTFREQ", "PMRA    ", "PMDEC   " 
       } ;

  char *tform[tfields] = 
     { "1J   ", "16A", "1J ", "4A ", "2E ", "2E ", "2E ", "2E ", "2D ", 
//       "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "2D ", "2D ", "1D ", "1D "
       "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D " // JNC JNC8apr16 CASA 4.5.2
     } ;

  char *tunit[tfields] =
     { "   ", "   ", "   ", "   ", "JY ", "JY ", "JY ", "JY ", 
       "HZ ", "HZ ", "DEGREES", "DEGREES", "YEARS  ", "DEGREES", 
       "DEGREES", "M/SEC  ", "HZ     ", "DEG/DAY", "DEG/DAY"
     } ;
  
  status=0;
  if ( ffcrhd(fptr,&status)) return status ;

  if (ffphbn (fptr, sources,tfields,ttype,tform,tunit,
              "AIPS SU",0,&status)) return status ;
  if (ffpkyj(fptr,"EXTVER", 1,"Version number of table",&status)) return status ;
  if (ffpkyj(fptr,"NO_IF" ,user->sidebands, "num_IF" ,&status)) return status ;
  if (ffpkys(fptr,"VELTYP" ,"TOPOCENT", "Velocity type" ,&status)) return status ;
  if (ffpkys(fptr,"VELDEF" ,"RADIO", " Velocity definition" ,&status))
    return status ;
  if (ffpkyj(fptr,"FREQID" ,1, "Frequency ID" ,&status)) return status ;

  for (row=1; row <= sources; row++)
  { 
    for (i=0; i<scans; i++) if (srec[i].source_id == row)break ;
    if (i == scans) 
    { fprintf(stderr,"FATAL: could not locate source %d\n", row) ;
      return -1 ;
    }
    source = &srec[i].scan->source ;
    sut.id = srec[i].source_id ;

    strcpy(sut.source, source->object) ;
    replace_nulls(sut.source, 16) ;
    sut.qual =  source->qual;
    for (j=0; j<4; j++)sut.calcode[j] = ' ';
    if(source->calcode>0)sut.calcode[0]='C';
    sut.bandwidth = user->channels*user->chan_inc*source->ch_width ;
    epoch=2000.0; //uvw coordinates now always rotated to 2000.0
    if(!user->recentre)
      app2j2000(source->ra_app,source->dec_app,user->recfile.mjd_ref,
		&ra,&dec);
    else
      app2j2000(burst->ra_app,burst->dec_app,user->recfile.mjd_ref,
		&ra,&dec);
    sut.ra   = ra*180/M_PI ;
    sut.dec  = dec*180/M_PI ;
    sut.epoch= epoch1;
    sut.ra_app = source->ra_app*180/M_PI ;
    sut.dec_app = source->dec_app*180/M_PI ;
    sut.pmra = (source->dra * 86400)*180/M_PI ;
    sut.pmdec = (source->ddec * 86400)*180/M_PI ;
    sut.lsrvel = source->lsrvel[0] ; //JNC8apr16 casa4.2.5
    sut.rest_freq= source->rest_freq[0] ; //JNC8apr16 casa4.2.5
    for (j=0; j<2; j++)  /* second entry is dummy */
    { //int pol=user->stokes_allow[0]/2;  /* 0=>130, 1=>175 */
      sut.iflux[j] = source->flux.i ; 
      sut.qflux[j] = source->flux.q ;
      sut.uflux[j] = source->flux.u ;
      sut.vflux[j] = source->flux.v ;
//      sut.lsrvel[j] = source->lsrvel[pol] ;
//      sut.rest_freq[j] = source->rest_freq[pol] ;
      sut.freq_off[j] = 0.0; /* All offsets in the FQ table */
    }
    if (corr->endian & LittleEndian)
    {
      swap_long(&sut.id, 1) ;
      swap_long(&sut.qual, 1) ;
      swap_long(sut.iflux, 4*2) ;
//      swap_d(sut.freq_off, 3*2+8) ;
      swap_d(sut.freq_off, 3*2+6) ; //JNC8apr16 casa4.2.5

    }
    fprintf(user->lfp,"Source %d : %12s\n",row,source->object) ;
    ffptbb(fptr, row, 1, sizeof(struct su_struct), (unsigned char *)&sut, &status) ;
  }
  return status ;
}
// create the header of the random group UV FITS file
int init_hdr( fitsfile *fptr, SvSelectionType *user, int pcount, int *status,
	      double epoch1)
{ int j ;
  int naxis, bitpix, simple,extend,blocked,groups;      
  long naxes[7],gcount=1 ;
  char ref_date[32], keynew[16] ;
  CorrType *corr = user->corr ;
  ProjectType *proj = &user->srec->scan->proj ;
  double epoch ;

  char *ctype[] = {"COMPLEX", "STOKES", "FREQ", "IF","RA","DEC"};
  float crval[] = {0.0,1.0,-1.0, 5.00E8,1.0,0.0,0.0};
  float cdelt[] = {0.0,1.0,-1.0, 4.8828125e4, 1.0,1.0E0,1.0E0};
  float crpix[] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  float crota[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  char *ptype[] = {"UU---SIN", "VV---SIN", "WW---SIN",
                     "BASELINE", "DATE", "DATE","SOURCE","FREQSEL"};
  float pscal[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  float pzero[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  int   gmrt2aips_stokes[]={-1.0,-3.0,-4.0,-2.0}; /* AIPS RL convention */
  int   gmrt2aips_stokes_xy[]={-5.0,-7.0,-8.0,-6.0}; /* AIPS XY convention */
  char  history[64];

  strncpy(ref_date , mjd2iau_date(corr->daspar.mjd_ref), 10) ;
  ref_date[10] = 0 ;

  if (pcount !=8)
    { fprintf (stderr, "INVALID PCOUNT = %d\n", pcount) ; return -1 ; }
  simple = TRUE;
  bitpix= -32;
  naxis = 7;
  naxes[0]=0;
  naxes[1]=3;
  naxes[2]= user->stokes ;
  naxes[3]= user->channels; //default: 1 channel per random group for SPOTLIGHT
  naxes[4]= user->sidebands ;
  naxes[5]=1;
  naxes[6]=1;
  epoch = 2000.0;// uvw now always in J2000
  crval[3] = RefFreq ;
  crpix[3] = 0.5 ;  /*GMRT frequency referes to the edge of the channel*/
  cdelt[3] = ChanWidth ;
  crval[2] = -1 ;    /* default RL since AIPS can't deal with linear */
  if(user->stokes_type==1) crval[2] = -5 ; /* user want's XY stokes labels */
  else crval[2] = -1 ;  /* assume RR LL RL LR order except for only 1 pol*/
  if(user->stokes==1)crval[2]=gmrt2aips_stokes[user->stokes_allow[0]];
  if(user->stokes==2)
  { int s0,s1;
    if(user->stokes_type==0)  /* RL */
    { s0=gmrt2aips_stokes[user->stokes_allow[0]];
      s1=gmrt2aips_stokes[user->stokes_allow[1]];
    }else    /* XY */
    {s0=gmrt2aips_stokes_xy[user->stokes_allow[0]];
     s1=gmrt2aips_stokes_xy[user->stokes_allow[1]];
    }
    crval[2]=s0;cdelt[2]=s1-s0;
  }
  extend= TRUE;
  blocked= TRUE;
  groups= TRUE;
  gcount = 1 ;
  proj->title[NAMELEN-1] = 0 ;
  if ( ffphpr (fptr,simple,bitpix,naxis,naxes,pcount, 
               gcount,extend,status) ) return *status ;
  if ( ffpkys(fptr,"OBJECT","MULTI", "",status) ) return *status ;
  if ( ffpkys(fptr,"TELESCOP","GMRT", "",status) ) return *status ;
  if ( ffpkys(fptr,"INSTRUME","GMRT", "",status) ) return *status ;
  if ( ffpkys(fptr,"OBSERVER",proj->observer, "",status) ) return *status ;
  //order is given as "TB" below, which is not correct at the moment - one
  //would need to do a uvsrt or equivalent to get the data in time-baseline
  //order
  if ( ffphis(fptr,"AIPS SORT ORDER = 'TB' ",status) ) return *status ;
  sprintf(history,"Created by SVFITS Version %5.3f",svVersion());
  if ( ffphis(fptr,history,status) ) return *status ;
  sprintf(history,"GMRT CORR Version %s",corr->version);
  if ( ffphis(fptr,history,status) ) return *status ;
  sprintf(history,"GMRT Project Code %s",proj->code);
  if ( ffphis(fptr,history,status) ) return *status ;
  if ( ffphis(fptr,history,status) ) return *status ;
  if ( ffpkys(fptr,"PROJECT",proj->title, "",status) ) return *status ;
  if ( ffpkys(fptr,"DATE-OBS",ref_date, "",status) ) return *status ;
  if ( ffpkys(fptr,"DATE-MAP",ref_date, "",status) ) return *status ;
  if ( ffpkyd(fptr,"BSCALE",1.0E+00 ,11,"",status) ) return *status ;
  if ( ffpkyd(fptr,"BZERO",0.0E+00 ,11,"",status) ) return *status ;
  if ( ffpkys(fptr,"BUNIT","UNCALIB", "",status) ) return *status ;
  if ( ffpkyd (fptr,"EPOCH",epoch, 9, "",status) ) return *status ;
  if ( ffpkye (fptr,"ALTRPIX",1.0, 9 , "",status) ) return *status ;
  if ( ffpkye (fptr,"ALTRVAL",0.0, 9 , "",status) ) return *status ; // Added ALTRVAL, value set to 0.0 JNC+SSK, 29th Apr 2021
  for (j=2; j < naxis+1; j++)
    {
      if ( ffkeyn ("CTYPE",j,keynew,status) ) return *status ;
      if (ffpkys(fptr,keynew,ctype[j-2], "",status)) return *status ;
      if ( ffkeyn ("CRVAL",j,keynew,status) ) return *status ;
      if (ffpkye(fptr,keynew,crval[j-1],10,"" ,status) ) return *status ;
      if ( ffkeyn ("CDELT",j,keynew,status) ) return *status ;
      if ( ffpkye (fptr,keynew,cdelt[j-1],9,"" ,status) ) return *status ;
      if ( ffkeyn ("CRPIX",j,keynew,status) ) return *status ;
      if ( ffpkye (fptr,keynew,crpix[j-1],9,"" ,status) ) return *status ;
      if ( ffkeyn ("CROTA",j,keynew,status) ) return *status ;
      if ( ffpkye (fptr,keynew,crota[j-1],10,"" ,status) ) return *status ;
    }
  for ( j=1; j <= pcount ; j++)
    {
      if (ffkeyn ("PTYPE",j,keynew,status) ) return *status ;
      if (ffpkys (fptr,keynew,ptype[j-1],"", status) ) return *status ;
      if (ffkeyn ("PSCAL", j,keynew, status) ) return *status ;
      if ( ffpkye(fptr,keynew,pscal[j-1],11,"",status)) return *status ;
      if (ffkeyn ("PZERO",j,keynew,status)) return *status ;
      if (ffpkye (fptr,keynew,pzero[j-1],11,"",status)) return *status ;
    }
  return *status ;
}

int printerror( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
  
  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
  if (!status) return status ;

  fprintf(stderr, "\n*** Error occurred during program execution ***\n");
  ffgerr(status, status_str);        /* get the error status description */
  fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);
  
  if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
      fprintf(stderr, "\nError message stack:\n");
      fprintf(stderr, " %s\n", errmsg);
      
      while ( ffgmsg(errmsg) )  /* get remaining messages */
	fprintf(stderr, " %s\n", errmsg);
    }
  
  return( status );
}
void set_coordtype(SvSelectionType *user)
{ CorrType *corr=user->corr;
  double mjd_ref=corr->daspar.mjd_ref;
  int i;

  for(i=0;i<strlen(user->coord_type);i++)
    user->coord_type[i]=tolower(user->coord_type[i]);
  /* user specification overides default */
  if(!strcmp(user->coord_type,"apparent"))
  { user->force_app=1; return ;}
  if(!strcmp(user->coord_type,"mean"))
  { user->force_app=0; return ;}
 
  /* apparent co-ordinates from 31/12/2000 to 15/08/2001 */
  if(mjd_ref > 51908.77 && mjd_ref < 52135.77) user->force_app=1;
  else user->force_app=0;

 return;

}
/*
  function to convert only the records in which there is a burst signal into
  random groups. If user->all_chan is set, then a multi-channel group is
  created and written to the FITS file (see make_allchan_group() for more
  details), otherwise, a single channel file is created, where only those
  channels which contain the burst signal are written to the FITS file.
  See copy_onechan_group() for more details.

  The user structure contains all the required meta data, fptr is a pointer
  to the output FITS file, which is assumed to be initiated, with the write
  pointer at the location where the groups are to be written.

  jnc may 2025
*/
int copy_burst(SvSelectionType *user,fitsfile *fptr){
  BpassType      *bpass=&user->bpass;
  RecFileParType *rfile=&user->recfile;
  int             rec_per_slice=rfile->rec_per_slice;
  int             baselines,channels;
  long            k,group_size,bufsize;
  int             recl,status=0;
  int             gcount,groups;
  int             b,i,start_rec,n_rec,rec,idx,slice,flagged;
  int             n_files,file_order[MaxRecFiles];
  char           *visbuf,*rbuf;
  
  // allocate space for raw visibilities (one full timeslice is processed
  // at a time). Space for the output random groups is allocated inside the
  // copy_* functions
  recl=user->corr->daspar.channels*user->corr->daspar.baselines*sizeof(float);
  bufsize=rec_per_slice*recl;
  if((rbuf=malloc(bufsize))==NULL){
    fprintf(stderr,"Malloc error for %ld bytes\n",bufsize);
    return -1;
  }
  // allocate memory for the mean spectrum and the amplitude bandpass
  channels=user->corr->daspar.channels; //input channels
  baselines=user->baselines;//output baselines
  for(b=0;b<baselines;b++){
    if((bpass->off_src[b]=(Complex*)malloc(channels*sizeof(Complex)))==NULL)
      {fprintf(stderr,"Malloc error in make_bpass\n"); return -1;}
    if((bpass->abp[b]=(float*)malloc(channels*sizeof(float)))==NULL)
      {fprintf(stderr,"Malloc error in make_bpass\n"); return -1;}
  }
  
  gcount=1;//FITS numbering starts from 1
  group_size=user->channels*user->stokes*sizeof(Cmplx3Type)+sizeof(UvwParType);

  // read visibilties from raw files, write out random groups to FITS FILE
  // loop over files and slices, and process one full slice at a time.
  if((n_files=get_file_order(user,file_order))<0)
    {fprintf(stderr,"No data found\n");return -1;}
  for(i=0;i<n_files;i++){
    idx=file_order[i];
    start_rec=rfile->start_rec[idx]; n_rec=rfile->n_rec[idx];
    //loop over slices, work with entire slice even if burst occupies only
    //part of the slice
    for(rec=start_rec;rec<start_rec+n_rec;rec+=rec_per_slice){
      slice=rec/rec_per_slice;
      if(!user->fake_data)
	if(read_slice(user,idx,slice,rbuf)) return -1;
      if((groups=copy_vis(user,idx,slice,start_rec,n_rec,rbuf,&visbuf))<0)
	{fprintf(stderr,"Error processing %s\n",rfile->fname[idx]);return -1;}
      if(!user->all_chan && user->do_flag){//flagging all chans is too expensive
	if((flagged=clip(visbuf,user,idx,slice,groups))<0)return -1;
	fprintf(user->lfp,"File: %d Slice:%d  Flagged %d of %d visibilities\n",
		idx,slice,flagged,groups);
      }
      k =  groups*group_size;
      if (user->corr->endian & LittleEndian)swap_long(visbuf, k/sizeof(float));
      ffptbb(fptr,gcount,1,(long)k,(unsigned char *)visbuf,&status) ; 
      if (status){printerror (status);break;}
      fprintf(stderr,"File %d Slice %d wrote %12.4e MBytes\n",idx,slice,
	      (1.0*k)/1.0e6);
      gcount += groups;
      free(visbuf); // allocated inside the copy_* functions
    }
  }
  gcount--; // total number of groups written

  for(b=0;b<baselines;b++){
    free(bpass->off_src[b]);
    free(bpass->abp[b]);
  }
  free(rbuf);
  
  return gcount;
}
/*
  function to copy all of the input raw visibilities into random groups,
  irrespective of whether the data contains a burst or not. This could be
  useful for making a continuum map for astrometry etc.

  The user structure contains all the required meta data, fptr is a pointer
  to the output FITS file, which is assumed to be initiated, with the write
  pointer at the location where the groups are to be written.

  jnc may 2025
 */
int copy_allvis(SvSelectionType *user,fitsfile *fptr){
  RecFileParType *rfile=&user->recfile;
  int             rec_per_slice=rfile->rec_per_slice;
  int             baselines=user->baselines;
  CorrType       *corr=user->corr;
  int             channels=corr->daspar.channels;
  int             stokes=user->stokes;
  long            k,group_size,bufsize;
  int             recl,gcount,groups,status=0;
  int             b,idx,slice,flagged;
  char           *visbuf,*rbuf;
  
  // allocate space for raw visibilities (one full timeslice is processed
  // at a time). Space for the output random groups is allocated inside the
  // copy_* functions
  recl=user->corr->daspar.channels*user->corr->daspar.baselines*sizeof(float);
  bufsize=rec_per_slice*recl;
  if((rbuf=malloc(bufsize))==NULL){
    fprintf(stderr,"Malloc error for %ld bytes\n",bufsize);
    return -1;
  }

  // allocate space for the output groups
  group_size=user->channels*user->stokes*sizeof(Cmplx3Type)+sizeof(UvwParType);
  groups=user->baselines/user->stokes;
  bufsize=groups*group_size;//only one record on output
  if((visbuf=(char*)malloc(bufsize))==NULL)
    {fprintf(stderr,"Malloc Error\n");return -1;}

  // read visibilties from raw files, write out random groups to FITS FILE
  // loop over files and slices, and process one full slice at a time.
  gcount=1;//FITS numbering starts from 1
  for(idx=0;idx<rfile->nfiles;idx++){
    for(slice=0;;slice++){//process data until we reach EoF
      if(read_slice(user,idx,slice,rbuf)<0) break; //Reached EoF
      if((groups=avg_vis(user,idx,slice,rbuf,visbuf))<0){
	fprintf(stderr,"Error processing %s slice %d\n",rfile->fname[idx],
		 slice);
	return -1;
      }
      k =  groups*group_size;
      if (user->corr->endian & LittleEndian)swap_long(visbuf, k/sizeof(float));
      ffptbb(fptr,gcount,1,(long)k,(unsigned char *)visbuf,&status) ; 
      if (status){printerror (status);break;}
      fprintf(stderr,"File %d Slice %d wrote %12.4e MBytes\n",idx,slice,
	      (1.0*k)/1.0e6);
      gcount += groups;
    }
  }
  gcount--; // total number of groups written

  free(rbuf);
  free(visbuf);

  return gcount;
}
int main(int argc, char **argv)
{ int             status=0,pcount=8,groups,flagged,idx,group_size;
  int             sources=1,frequencies=1; // for SU and FQ tables
  float           version;
  long            gcount;
  fitsfile       *fptr;
  char            outfile[LINELEN] ;
  SvSelectionType user;
  InitHdrType    *hdr;
  SourceParType  *source;
  char            antfile[PATHLEN],uparfile[PATHLEN];
  int             c;
  extern char    *optarg;
  extern int      optind,opterr;

  strcpy(antfile,"antsamp.hdr");
  strcpy(uparfile,"svfits_par.txt");
  if(argc-1){
    while((c=getopt(argc,argv,"a:u:h"))!=-1){
      switch(c){
	case 'a': strncpy(antfile,optarg,PATHLEN-1); break;
	case 'u': strncpy(uparfile,optarg,PATHLEN-1);break;
        case 'h': // default to next case
        case '?': fprintf(stderr,
		  "Usage: svfits [-a AntSampFile] [-u UserParmFile]\n");
        	  return -1;
      }
    }
  }
 
  fprintf(stdout,"     ----- SVFITS  version %.2f  -----  \n",svVersion()) ;
  fprintf(stdout,"FITSIO version number = %f\n",ffvers(&version));

  //need to initialize user, corr, etc over here
  user.hdr=(InitHdrType*)malloc(sizeof(InitHdrType));
  hdr=user.hdr;
  hdr->scans=1;
  user.srec = (ScanRecType *) malloc(hdr->scans*sizeof(ScanRecType)) ;
  user.srec->scan=(ScanInfoType *) malloc(hdr->scans*sizeof(ScanInfoType)) ;
  source=&user.srec->scan->source;
  user.corr=(CorrType*)malloc(sizeof(CorrType));
  if(init_user(&user,uparfile,antfile)<0){return -1;}
  RefFreq   = source->freq[0];
  ChanWidth = source->ch_width*source->net_sign[0] ;

  //Initialize the FITS file and write the main header
  strcpy(outfile,user.fitsfile);
  { char *p ; // remove file if it already exists
    if ((p = rindex(outfile,'/')) == NULL)p = outfile ;
    while (*p){*p = toupper(*p) ; p++ ; }
    remove(outfile) ;
  }
  if ( ffinit(&fptr,outfile,&status)) return printerror(status);
  if (init_hdr(fptr,&user,pcount,&status,user.epoch) != 0)
  { if (status) return printerror(status);
    fprintf(stderr,"FATAL Error in init_hdr\n");
    return -1 ;
  }

  if(user.all_data)
    {if((gcount=copy_allvis(&user,fptr))<0) return -1;}
  else
    {if((gcount=copy_burst(&user,fptr))<0) return -1;}

  // update total number of groups   
  if(ffukyj(fptr,"GCOUNT",gcount," ",&status)) return printerror(status);
  if(ffclos(fptr,&status)) return printerror(status) ;

  // create the FITS file tables
  if ( ffopen (&fptr, outfile, READWRITE,&status)) printerror(status) ;
  if (writeagtable(fptr,&user,user.iatutc))//antenna table
    printerror(status);
  if (writefqtable(fptr,&user,user.scans,frequencies))//frequency table
    printerror(status) ;
  if (writesutable(fptr, &user, user.scans, sources,user.epoch))//source table
    printerror(status) ;
  if(ffclos(fptr,&status)) printerror( status );
  
  free(user.hdr) ;
  fprintf(stdout,"Created %s\n",user.fitsfile);
  fprintf(stdout,"See svfits.log for details\n");

  return 0 ;
}


