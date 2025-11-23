#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<time.h>
#include<omp.h>
#include<errno.h>

#include"fitsio.h"
#include "gmrt_newcorr.h"
#include "svio.h"


enum{Hdrlen=81920,Len=80};
static char Buf[Hdrlen], *Bufp ;// local buffer for reading hdr files
const char cont_char = '*', begin_char = '{', end_char = '}';

// The enum and the char array need to match exactly, and are used to
// parse the input parameter file.
enum {NFiles,Path,InputFile,HaveIdx,Lta,NLta,FitsFile,AntMask,AllChan,AllData,
      NCHAV,IATUTC,Epoch,ObsMjd,FreqSet,RaApp,DecApp,RaMean,DecMean,
      StokesType,BurstName,BurstMJD,BurstTime,BurstDT,BurstIntWd,BurstWd,
      BurstDM,BurstDDM,BurstFreq,BurstBmId,BurstRa,BurstDec,UpdateBurst,
      DoFlag,Thresh,DoBand,DoBase,PostCorr,DropCsq,NumThreads,Recentre,
      CoordType,SvVars} SvVarType ;
char *SvVarname[SvVars] =
  {"NFILE","PATH","INPUT","HAVE_IDX","LTA","N_LTA","FITS","ANTMASK","ALL_CHAN",
   "ALL_DATA","NCHAV","IATUTC","EPOCH","OBS_MJD","FREQ_SET",
   "RA_APP","DEC_APP","RA_MEAN","DEC_MEAN","STOKES_TYPE","BURST_NAME",
   "BURST_MJD","BURST_TIME","BURST_DT","BURST_INTWD","BURST_WIDTH","BURST_DM",
   "BURST_DDM","BURST_FREQ","BURST_BM_ID","BURST_RA","BURST_DEC","UPDATE_BURST",
   "DO_FLAG","THRESH","DO_BAND","DO_BASE","POST_CORR","DROP_CSQ",
   "NUM_THREADS","RECENTRE","COORD_TYPE"};
static double K0=4.15e-3; // DM constant in seconds

float svVersion(void){//sets the version number for log and history
  return 0.972;
}

/*
  parse the inputs supplied by the user in the parameter file specified by the
  -u option. These inputs will override the defaults set in init_user(). The
  svio.h header file and the template user input file has some information on
  what the different parameters are mean for.

  jnc apr 25
*/
int svuserInp (char *filename, SvSelectionType *user ){
  CorrType      *corr=user->corr;
  DasParType    *daspar=&user->corr->daspar;
  SourceParType *source=&user->srec->scan->source;
  char           str[2048],key_str[32] ;
  int            val[2048];
  FILE          *fp;
  int            var, np, k ;

  if((fp= fopen(filename, "rt")) ==NULL)
    {fprintf(stderr,"Cannot open %s\n",filename); return -1;}
  
  while ( fgets(str, 2040, fp) )
  { char *p = str ;
    while (*p && (*p == ' '))p++ ;
    if (*p == 0)continue ;
    if ( (*p == '#') || (*p == '*') )continue ;
    key_str[0] = 0 ;
    sscanf(p, "%s",key_str) ;
    if(!strlen(key_str)) continue;
    for (var=0; var < SvVars; var++)
      if (!strncasecmp(key_str,SvVarname[var],strlen(SvVarname[var])))
	  break ;
    if (var == SvVars)
      {fprintf(stderr,"Unrecognized Parameters %s\n",key_str); return -1;}

    while (*p && !isspace((int)*p))p++ ;
    while (*p && isspace((int)*p))p++ ;
    for(k=strlen(p);k>=0;k--)if(p[k]=='!'){p[k]='\0';break;}
    switch(var){
      int i,nchan;
      double freq0,freq1;
      char *f;
      case NFiles   :sscanf(p,"%d",&user->recfile.nfiles); break;
      case Path     :strncpy(user->recfile.path,p,PATHLEN-1); break;
      case InputFile:f=strtok(p,",");i=0;
                      while(f!=NULL && i < user->recfile.nfiles){
	 	      strcpy(user->recfile.fname[i],f);f=strtok(NULL,",");
		      i++;
		     }
		     break;
      case HaveIdx  :sscanf(p,"%d",&user->recfile.have_idx);break;
      case Lta      :sscanf(p,"%lf",&user->lta);break;
      case NLta   :sscanf(p,"%d",&user->n_lta);break;
      case AntMask  :sscanf(p,"%u",&user->antmask);break;
      case AllChan  :sscanf(p,"%d",&user->all_chan);break;
      case AllData  :sscanf(p,"%d",&user->all_data);break;
      case NCHAV    :sscanf(p,"%d",&user->nchav);break;
      case FitsFile :sscanf(p,"%s", user->fitsfile);break ;
      case IATUTC   :sscanf(p,"%f",&user->iatutc);break ;
      case ObsMjd   :sscanf(p,"%lf",&user->corr->daspar.mjd_ref);break ;
      case FreqSet  :i=sscanf(p,"%lf:%lf:%d\n",&freq0,&freq1,&nchan);
	             if(i!=3)
		       {fprintf(stderr,"Format error in FREQ_SET %d\n",i); return -1;}
        	     source->freq[0]=source->freq[1]=freq0; 
                     source->ch_width=fabs(freq1-freq0)/nchan;
		     if(freq1<freq0)
		       source->net_sign[0]=source->net_sign[1]=-1;
		     else
		       source->net_sign[0]=source->net_sign[1]=1;
                     corr->daspar.channels=nchan;
		     break;
      case Epoch    :sscanf(p,"%lf",&user->epoch);break;
      case RaApp    :sscanf(p,"%lf",&user->srec->scan->source.ra_app);break;
      case DecApp   :sscanf(p,"%lf",&user->srec->scan->source.dec_app);break;
      case RaMean   :sscanf(p,"%lf",&user->srec->scan->source.ra_mean);break;
      case DecMean  :sscanf(p,"%lf",&user->srec->scan->source.dec_mean);break;	      
      case StokesType: sscanf(p,"%hd",&user->stokes_type);break;
      case BurstName:strcpy(user->burst.name,p);
	             strcpy(user->srec->scan->source.object,user->burst.name);
		     break;
      case BurstMJD :sscanf(p,"%lf",&user->burst.mjd);break;
      case BurstTime:sscanf(p,"%lf",&user->burst.t);break;
      case BurstDT  :sscanf(p,"%lf",&user->burst.dt);break;
      case BurstIntWd:sscanf(p,"%lf",&user->burst.int_wd);break;
      case BurstWd  :sscanf(p,"%lf",&user->burst.width);break;
      case BurstDM  :sscanf(p,"%lf",&user->burst.DM);break;
      case BurstDDM :sscanf(p,"%lf",&user->burst.dDM);break;
      case BurstFreq:sscanf(p,"%lf",&user->burst.f);break;
      case BurstBmId:sscanf(p,"%d",&user->burst.bm_id);break;
      case BurstRa  :sscanf(p,"%lf",&user->burst.ra_app);break;
      case BurstDec :sscanf(p,"%lf",&user->burst.dec_app);break;
      case UpdateBurst:sscanf(p,"%d",&user->update_burst);break;
      case Recentre :sscanf(p,"%d",&user->recentre);break;
      case DoFlag   :sscanf(p,"%d",&user->do_flag);break;
      case Thresh   :sscanf(p,"%f",&user->thresh);break;
      case DoBand   :sscanf(p,"%d",&user->do_band);break;
      case DoBase   :sscanf(p,"%d",&user->do_base);break;
      case PostCorr :sscanf(p,"%d",&user->postcorr);break;
      case DropCsq  :sscanf(p,"%d",&user->drop_csq);break;
      case CoordType:break; // depreciated
      case NumThreads:sscanf(p,"%d",&user->num_threads);break;
    }
  }

  fclose(fp);
  return 0;
}
/* set up some booking keeping information for each baseline needed at the
   time of copying the data, conversion to FITS etc. Sets up the map between
   the input baseline order and the output random group order. The map is
   setup so that baselines between the same antenna pair, but for different
   stokes parameters will become adjacent (as required) in the output random
   group.

   jnc mar 2025
   28mar25 fixed bug in loop order

   corrected the offset to reflect that there is no timestamp at the start
   of every record. Instead there is a timestamp at the start of every
   rec_per_slice records. The timestamp is read separately and so does not
   need to be reflected in the offset here.
*/
int init_vispar(SvSelectionType *user){

  InitHdrType   *hdr=user->hdr;
  ScanInfoType  *scan=user->srec->scan;
  CorrType      *corr=user->corr;
  VisParType    *vispar=&user->vispar;
  VisInfoType   *visinfo=vispar->visinfo;
  int            antennas=MAX_ANTS;
  SamplerType   *samp=corr->sampler;
  SourceParType *source=&scan->source;
  unsigned int   antmask=user->antmask; // antennas whose data to copy
  int            baselines=corr->daspar.baselines;
  BaseParType   *base=corr->baseline;
  unsigned int   stride=corr->daspar.channels*sizeof(float); //half-float re,im
  unsigned long  off;
  int            ant0,ant1,band,bs,bs1,recl;

  recl=corr->daspar.channels*sizeof(float); //datasize for one baseline
  vispar->vis_sets=user->stokes;

  // organize so that the output visibilities are in antenna collation order
  // and that baselines with the same antenna pair but different stokes are
  // next to each other.  The sign of the imaginary component is set depending
  // on the antenna index order as well as whether it is USB/LSB
  for(bs1=0,ant0=0;ant0<antennas;ant0++){
    if(!(1<<ant0&antmask)) continue;
    for(ant1=ant0+1;ant1<antennas;ant1++){//avoid self
      if(!(1<<ant1&antmask)) continue;
      for(band=0;band<user->stokes;band++){
	for(bs=0;bs<baselines;bs++){
	  int a0,a1,b0,b1;
	  a0=base[bs].samp[0].ant_id;b0=base[bs].samp[0].band;
	  a1=base[bs].samp[1].ant_id;b1=base[bs].samp[1].band;
	  if((b0!=band) ||(b1!=band)) continue;
	  if(ant0==a0 && ant1==a1){
	    if(source->net_sign[0]<0) visinfo[bs1].flip=0; 
	    else  visinfo[bs1].flip=1; 
	    visinfo[bs1].ant0=ant0;visinfo[bs1].ant1=ant1;
	    visinfo[bs1].band0=b0;visinfo[bs1].band1=b1;
	    visinfo[bs1].drop=0; 
	    visinfo[bs1].off=bs*recl;
	    bs1++;break;
	  }else{
	    if(ant0==a1 && ant1==a0){
	      if(source->net_sign[0]<0) visinfo[bs1].flip=1;
	      else visinfo[bs1].flip=0;
	      visinfo[bs1].ant0=ant0;visinfo[bs1].ant1=ant1;
	      visinfo[bs1].band0=b0;visinfo[bs1].band1=b1;
	      visinfo[bs1].drop=0; 
	      visinfo[bs1].off=bs*recl;
	      bs1++;break;
	    }
	  }
	}
      }
    }
  }

  return bs1;
}
/*
  various string handling functions used by the other functions in this file.
  taken directly from gvfits and/or csubs1.c
*/
int str_cmp(const char *s1, const char *s2)
{
  while (*s1 && toupper(*s1) == toupper(*s2)) { s1++; s2++; }
  return *s1 - *s2;
}
int search_str(char *s1, char **s2, int n)
{ int i;
  for (i=0; i < n; i++)
    if (str_cmp(s1, s2[i]) == 0) break;
  return i;
}
char *next_val(char *s)
{ static char *p;
  char *q;
  if (s != NULL) p = s;
  while (isspace(*p)) p++;
  if (*p == 0) return NULL;
  q = p;
  if (*p == ':') p++; else while (isdigit(*p)) p++;
  return q;
}
int get_range(char *str, short *arr, int max)
{ int i, j, k, l, m, n=0, ind=0, val[100] ;
  char *p ;

  p = next_val(str);
  while (p != NULL)
  { if (sscanf(p,"%d",&val[n++]) != 1)
    { fprintf(stderr,"get_range: Illegal range\n"); n--; break; }
    p = next_val(NULL); if (p == NULL) break;
    if (*p == ':') { val[n++] = INT_MIN ; p = next_val(NULL) ; }
  }
  if (val[n-1] == INT_MIN) n--;
  if (n <= 0) return 0;
  val[n] = 0;
  for (i=0, m=0; i < n; i = ind+1)
  { k = 1;
    ind = i ;      /* fixed on Sep 9, 1995 ! */
    if (i+2 < n && val[i+1] == INT_MIN)
    { j = val[ind = i+2];
      if (i+4 < n && val[i+3] == INT_MIN) k = val[ind = i+4];
    }
    else j = val[i];
    if (j >= max) j = max-1;
    for (l=val[i]; l <= j; l += k) arr[m++] = l;
  }
  return m;
}
void copy_text(FILE *f)
{ char *p;
  while (fgets(Bufp, Len, f))
  { p = Bufp;
    while (*Bufp) Bufp++;
    if (*p == end_char) break;
  }
}
/*
  get the antenna configuration, i.e. names, positions etc. This function
  reads these parameters from an input text file (the default name is
  hardcoded, and is in the same format as used by all the other programs in
  the das pipeline).  For realtime pipelined operations however, one should
  take the antenna table from the "corr" structure which (along with the
  "scan" structure) should be dumpled by the raw visibility writing program.
  
*/
void get_antenna(FILE *f, AntennaParType *antenna)
{ char           *p;
  int             seq, j, id;
  double          fixed, bx,by,bz;
  AntennaParType *ant;

  for (seq=0; seq < MAX_ANTS; seq++){
    for (j=0; j < MAX_BANDS; j++){
      antenna[seq].samp_id[j] = MAX_SAMPS;// reverse look filled in get_sampler
      antenna[seq].p0[j] = 0; // per antenna phase set to zero
    }
  }

  seq = 0;
  while (fgets(Bufp, Len, f))
  { p = Bufp;
    while (*Bufp) Bufp++;
    if (*p == cont_char) continue;
    if (*p == end_char) break;
    if (strncmp(p, "ANT", 3) != 0) continue;
    sscanf(p+3,"%d", &id); if (id >= MAX_ANTS) continue;
    ant = &antenna[id];
    sscanf(p+10,"%s %lf %lf %lf %lf", ant->name, &bx, &by, &bz, &fixed);
    ant->bx = bx ;      ant->by = by ;    ant->bz = bz ;
    if (id < seq) fprintf(stderr,"Antenna Sequence changed..\n");
    seq = id+1;
    for (j=0; j < MAX_BANDS; j++) ant->d0[j] = fixed;
  }
}

void get_bands(FILE *f, char **band_name)
{ int i=0, id;
  char *p;
  while (fgets(Bufp, Len, f))
  { p = Bufp;
    while (*Bufp) Bufp++;
    if (*p == cont_char) continue;
    if (*p == end_char) break;
    if (strncmp(p, "BAND", 4) != 0) continue;
    sscanf(p+4,"%d", &id);
    sscanf(p+10,"%s", band_name[id]);
    if (id != i) fprintf(stderr,"Band Sequence changed..\n");
    i++;
  }
}
/* get sampler connections from header file. Bare bones version for SPOTLIGHT
 raw data FITS conversion. This function reads these parameters from an input
 text file (the default name is hardcoded, and is in the same format as used
 by all the other programs in the das pipeline).  For realtime pipelined
 operations however, one should take the antenna table from the "corr"
 structure which (along with the "scan" structure) should be dumpled by
 the raw visibility writing program.

*/
int get_sampler(FILE *f, CorrType *corr){
  int             ant_id, seq=0, id, band,i,j;
  char            str1[16],str2[16], *p;
  AntennaParType *antenna = corr->antenna;
  SamplerType    *samp    = corr->sampler;
  DasParType     *daspar  =&corr->daspar ;
  char           *band_name[MAX_BANDS]; // copy for search_str

  band_name[0] =  corr->bandname[0];
  band_name[1] =  corr->bandname[1];
  band_name[2] =  corr->bandname[2];
  band_name[3] =  corr->bandname[3];

  daspar->samplers = 0 ;
  daspar->bandmask = 0 ;
  daspar->antmask  = 0 ;
  for (id=0; id < MAX_SAMPS; id++){ 
    samp[id].ant_id = MAX_ANTS;
    samp[id].band   = MAX_BANDS ;
    samp[id].fft_id = MAX_FFTS;
    samp[id].dpc    = MAX_SAMPS ;
  }

  seq = 0 ;
  while ((p = fgets(Bufp, Len, f)) != NULL){
    Bufp += strlen(Bufp);
    if (*p == cont_char) continue;
    if (*p == end_char) break;

    if (strncmp(p, "SMP", 3) != 0) continue;
    sscanf(p+3,"%d", &id);
    if (id >= MAX_SAMPS){fprintf(stderr,"Band SampId %d\n",id); continue;}
    sscanf(p+10,"%s %s", str1, str2 );
    /* antenna_name, band_name  to which it is connected */
    for (ant_id=0; ant_id < MAX_ANTS; ant_id++)
      if (str_cmp(str1, antenna[ant_id].name) == 0) break ;
    if (ant_id == MAX_ANTS ) continue ;
    band = search_str(str2, band_name, MAX_BANDS);
    if (band == MAX_BANDS){fprintf(stderr,"Bad Band %s\n",str2);continue;}
    daspar->antmask  |= (1 << ant_id);
    daspar->bandmask |= (1 << band);
    samp = corr->sampler + daspar->samplers ;
    samp->ant_id = ant_id;
    samp->band=band; // JNC 18Mar25 (?not there earlier??)
    antenna[ant_id].samp_id[band]=daspar->samplers ; /*reverse lookup */
    daspar->samplers++;
  }
  return 0;
}
/*
  stripped down version to fill up the antennas and bands (i.e. polarizations)
  corresponding to a given baseline. The order here is as per the clarification
  provided by Harshavardhan Reddy on 20/Mar/25 for the structure of SPOTLIGHT
  raw visibility dumpfiles.

  jnc mar 2025
*/
int get_baseline(SvSelectionType *user){
  CorrType     *corr=user->corr;
  unsigned int  antmask=corr->daspar.antmask;
  SamplerType  *samp=corr->sampler;
  BaseParType  *base=corr->baseline;
  int           s0,s1,b,band;

  if(user->stokes!=2)// raw dump is for stokes RR,LL only
    {fprintf(stderr,"Illegal Nstokes %d [Legal 2]\n",user->stokes); return -1;}

  // we have all RR baselines first, followed by all LL baselines
  for(b=0,band=0;band<user->stokes;band++){
    for(s0=0;s0<corr->daspar.samplers;s0++){
      int ant0=samp[s0].ant_id,band0=samp[s0].band;
      if(band0 !=band || !(1<<ant0 & antmask)) continue;
      for(s1=s0;s1<corr->daspar.samplers;s1++){
	int ant1=samp[s1].ant_id,band1=samp[s1].band;
	if(band1!=band || !(1<<ant1 & antmask)) continue;
	memcpy(&base[b].samp[0],samp+s0,sizeof(SamplerType));
	memcpy(&base[b].samp[1],samp+s1,sizeof(SamplerType));
	b++;
      }
    }
  }
  corr->daspar.baselines=b;

  return 0;
}

/*
  initialize the corr structure which carries informatio both about the
  correlator capacities, the settings chosen by the user, the antenna to
  sampler mapping, etc. This meta data in this structure is used in all
  instances when information about the visibilties is needed.
  
 */
int init_corr(SvSelectionType *user, char *hdrfile){
  char *p, str[32];
  FILE *f;
  int i ;
  CorrType    *corr=user->corr;
  CorrParType *corrpar=&corr->corrpar;
  DasParType  *daspar=&corr->daspar;

  //Historical USB-130 is RR, USB-175 is LL, LSB is unused
  char *band_name[MAX_BANDS] = 
    { "USB-130", "USB-175", "LSB-130", "LSB-175" } ;
  //2 stokes, 4 stokes (full polar) or 1 stokes (single pol)
  char *mac_modes[MACMODES] = {"RRLL","RRRL","RR" } ;

  //default correlator parameters. Only those relevant to the SPOTLIGHT
  //correlator are set.
  bzero(corr,sizeof(CorrType));
  strcpy(corr->version, "CORR-SPOTLIGHT-RAW");corr->version[NAMELEN-1]='\0';
  corr->endian=LittleEndian;// x86 is little endian
  for (i=0; i<MAX_BANDS; i++) strcpy(corr->bandname[i],band_name[i]) ;
  corrpar->pols=2;        //2 stokes
  corrpar->macmode=RRLL;//stokes are "RR" and "LL" (RRLL is an enum)
  corrpar->clock=400e6;   // default 200 MHz bandwidth
  corrpar->clksel=0; // unused here, in principle clock_used =clock/2<<clk_sel
  corrpar->channels=4096; 
  corrpar->sta=1; // adhoc  does not make a different to FITS conversion
  // duration of 1 sta in sec(changed for SPOTLIGHT)
  user->statime=2.0*corrpar->channels*corrpar->sta/corrpar->clock;
  corrpar->f_step=corrpar->clock/(corrpar->channels*2);// +ve 

  // default data acquisition parameters - (in principle not all correlator
  // output need be acquired, user could select subset).
    // antmask==1<<ant_id ==> include antenna index ant_id in antenna file
  daspar->antmask=4294967295; // 32 antennas (updated in get_sampler)
  daspar->bandmask=3; // RR, LL polarizations 
  daspar->samplers=64; // 64 samplers (updated in get_sampler) 
  daspar->baselines=MAX_BASE;// assumes 4 stokes (updated in get_baseline)
  daspar->channels=corrpar->channels;
  for(i=0;i<MAX_CHANS;i++)daspar->chan_num[i]=i;
  daspar->lta=64;   //sta cycles per lta == 1.30172 ms
  daspar->t_unit=1; //backwards compatibility
  daspar->mjd_ref=60676;// 1Jan25 (place holder, will be overwritten)
  daspar->dut1=0.0; //should be overwritten by BulletinA vals (but unused)
  
  // get the sampler and antenna parameters from the header file
  Bufp = Buf ;
  f = fopen(hdrfile, "rt");
  if (f == NULL){fprintf(stderr,"cannot open %s\n",hdrfile); return -1;}
  while (fgets(Bufp, Len, f))
  { if (*Bufp == cont_char) continue;
    while (isspace(*Bufp)) Bufp++;
    if (*Bufp == 0) continue;
    p = Bufp; Bufp += strlen(Bufp);
    if (*p == begin_char)
    { p++; sscanf(p,"%s",str);
      if (str_cmp(str,"Antenna.def") == 0)
        get_antenna(f, corr->antenna);
      else if (str_cmp(str,"Sampler.def") == 0)
      {get_sampler(f, corr);}
      else copy_text(f); //skip over textt
    }
    else if (strncmp(p,"END_OF",6) == 0) break;
  }
  fclose(f);

  if(get_baseline(user)) return -1;
  return 0;
}

/*
  function to compute the burst time and burst width, given the burst->mjd,
  the intrinsic width and the DM. The IST date and time of the burst is
  returned in burst->date, the IST time is returned in burst->t (seconds)
  and the burst width is returned in burst->width (sec). Typically in pipeline
  mode the burst->mjd and the burst->int_wd are all that are available.

   JNC Mar2025

*/
int update_burst(SvSelectionType *user){
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq=source->freq[0]/1.0e9; //GHz
  double           cw=source->ch_width/1.0e9;//GHz (positive)
  double           fb=burst->f/1.0e9;//GHz
  double           int_wd=burst->int_wd,DM=burst->DM;
  double           mjd=burst->mjd,f_mjd,i_mjd;
  CorrType        *corr=user->corr;
  double           fh,fl,th,tl,f_sec,tsec;
  int              hh,mm,ss;
  time_t           unix_time;

  if(source->net_sign[0]>0){fh=freq+corr->daspar.channels*cw;fl=freq;}
  else{fh=freq; fl=freq-corr->daspar.channels*cw;}

  //total observed width (including intrinsic width and dispersion)
  th=-int_wd/2+K0*DM*(1.0/(fh*fh)-1.0/(fb*fb));//burst start at highest freq
  tl= int_wd/2+K0*DM*(1.0/(fl*fl)-1.0/(fb*fb));//burst end at lowest freq
  burst->width=tl-th; // seconds

  //convert mjd to IST
  i_mjd=floor(mjd);f_mjd=mjd-i_mjd;
  //integer seconds since midnight 01/01/1970
  unix_time=(i_mjd-40587.0)*86400.0+floor(f_mjd*86400); 
  tsec=f_mjd*86400;
  f_sec=tsec-floor(tsec);//fractional seconds
  ctime_r(&unix_time,burst->date);//in IST now
  burst->date[strcspn(burst->date,"\n")]='\0';//remove trailing newline,if any
  sscanf(burst->date,"%*s %*s %*s %d:%d:%d %*s",&hh,&mm,&ss);
  burst->t=hh*3600+mm*60+ss+f_sec; //IST seconds
  if(user->do_log){
    fprintf(user->lfp,"Burst Date: %s\n",burst->date);
    fprintf(user->lfp,"Burst Time: %12.6f [%d %d %10.7f]\n",burst->t,hh,mm,
	    ss+f_sec);
    fprintf(user->lfp,"Burst Wdth: %12.6f\n",burst->width);
  }
  return 0;
}
/*
  initialize the parameters in the "user" structure with sensible values.
  This is the fundamental structure passed to all functions, and carries all
  of the meta data needed for processing the visibilities. This includes both
  the correlator and observation settings, the antenna to sampler mapping,
  antenna positions, as well as the selections made by the user for processing
  the visibilities into a FITS file.

  Many of the parameters (particularly those which the user specifies for the
  data format conversion from raw visibility dump to FITS) will be
  overwritten by the actual parameters used at run time. However it is helpful
  to initialize to sensible values since this allows for the code to be
  run and debugged without having any input files.

  jnc march 2025
*/
int init_user(SvSelectionType *user, char *uparfile, char *antfile,
	      char *bhdrfile, char *bulletinA){
  fprintf(stderr,"DEBUG init_user: start, user=%p\n", (void*)user); fflush(stderr);

  fprintf(stderr,"DEBUG init_user: user->corr=%p\n", (void*)user->corr); fflush(stderr);
  VisParType       *vispar=&user->vispar;
  CorrType         *corr=user->corr;
  fprintf(stderr,"DEBUG init_user: corr=%p\n", (void*)corr); fflush(stderr);

  CorrParType      *corrpar=&corr->corrpar;
  DasParType       *daspar=&corr->daspar;
  RecFileParType   *rfile=&user->recfile;
  fprintf(stderr,"DEBUG init_user: user->srec=%p\n", (void*)user->srec); fflush(stderr);

  ScanRecType      *srec=user->srec;
  fprintf(stderr,"DEBUG init_user: srec=%p, srec->scan=%p\n", (void*)srec, (void*)srec->scan); fflush(stderr);

  ScanInfoType     *scan=srec->scan;
  fprintf(stderr,"DEBUG init_user: scan=%p\n", (void*)scan); fflush(stderr);

  SourceParType    *source=&scan->source;
  BurstParType     *burst=&user->burst;
  int               i;

  fprintf(stderr,"DEBUG init_user: pointers OK, opening log\n"); fflush(stderr);

  //open the log file
  if((user->lfp=fopen("svfits.log","w"))==NULL)
  { fprintf(stderr,"Cannot open svfits.log\n"); return -1;}

  //set timestamp to the middle of the integration
  user->timestamp_off=0.5*daspar->lta*user->statime;
  user->iatutc = 37.0; // default, leap seconds are probably now fixed
  user->epoch = 2000.0;// default is J2000 coordinates
  user->stokes=2; // RR, LL by default
  /* Set default stokes for RR,LL */
  user->stokes_allow[0] = 0 ;   // RR
  user->stokes_allow[1] = 3 ;   // LL
  user->stokes_allow[2] = user->stokes_allow[3] = -1 ;
  user->stokes_type=0; /* FITS stokes labels are RL */
  user->scans       = 1;
  user->sidebands   = 1;
  user->all_chan    = 0;// copy all channels for records with burst
  user->all_data    = 0;// copy all data, ignore burstpar
  user->nchav       = 1;// channels to average together (for all_data=1 only)
  user->start_chan  = 0;
  user->chan_inc    = 1;
  user->force_app   = 1; // force coordinates to be treated as apparent
  user->fake_data   = 0; // FOR DEBUGGING!!!! SET TO 0 FOR ACTUAL USE
  user->do_log=20;
  user->do_flag=1;// flag visibilities with outlier amplitudes
  user->update_burst=1;//recompute parms with some as input (see update_burst())
  user->thresh=5.0;// threshold for MAD flagging
  user->do_band=0; // do not do amplitude bandpass
  user->do_band=0; // do not subtract baseline (i.e. mean offsrc visibility)
  user->postcorr= 0;//generate a postcorrelation array output (debug use)
  user->drop_csq= 1;//drop csq baselines while making postcorr beam
  user->num_threads=16;//number of threads to use in parallel sections
  user->bpass.file_idx=-1;
  user->bpass.slice=-1;
  strcpy(user->fitsfile,"TEST.FITS") ;
  /* default parameters, reset by reading input scanfile*/
  fprintf(stderr,"DEBUG init_user: user->hdr=%p\n", (void*)user->hdr); fflush(stderr);
  fprintf(stderr,"DEBUG init_user: about to write user->hdr->scans\n"); fflush(stderr);
  user->hdr->scans=user->scans;
  fprintf(stderr,"DEBUG init_user: hdr->scans done, opening log again\n"); fflush(stderr);
  if((user->lfp=fopen("svfits.log","w"))==NULL)
    {fprintf(stderr,"Unable to open svfits.log\n"); return -1;}
  fprintf(stderr,"DEBUG init_user: calling init_corr with antfile=%s\n", antfile); fflush(stderr);
  // initialize the correlator settings
  init_corr(user,antfile); // hardcoded for now
  fprintf(stderr,"DEBUG init_user: init_corr done\n"); fflush(stderr);
  user->channels=1; //only one output channel by default
  user->antmask=1073741823;//30 antennas (C07 and S05 dropped)
  srec->corr=user->corr;
  srec->scannum=0;
  srec->scan_id=1;
  srec->source_id=1;
  srec->freq_id=1;
  strcpy(scan->proj.code,"TEST"); // replace with GTAC project code
  strcpy(scan->proj.observer,"DUMMY"); // replace with GTAC observer
  strcpy(scan->proj.title,"SPOTLIGHT"); // replace with GTAC observer
  scan->proj.antmask=daspar->antmask;
  scan->proj.bandmask=daspar->bandmask;
  scan->proj.seq=0;
  source->antmask=daspar->antmask;
  source->bandmask=daspar->bandmask;
  source->ch_width=corrpar->f_step;
  source->freq[0]=scan->source.freq[1]=500e6; // 500 MHz
  source->net_sign[0]=scan->source.net_sign[1]=1;//frq increases with chan
  source->ra_app=source->ra_mean=M_PI/6.0;
  source->dec_app=source->dec_mean=M_PI/4.0;
  strcpy(scan->source.object,"SPOTLIGHT");
  //set up some default burst parameters
  strcpy(burst->name,"TEST-BURST");
  burst->t=100.0;// 100 seconds after start
  burst->dt=0.0;
  burst->mjd=corr->daspar.mjd_ref + burst->t/86400;
  burst->DM=30.0;// change burst->width also if this is changed
  burst->dDM=0;  // error in DM (unused at the moment)
  burst->f=700.0e6; // Hz
  burst->width=0.23; //observed burst width (changes with burst->DM)
  burst->int_wd=1.0e-3;//intrinsic burst width
  burst->bm_id=0;
  burst->ra_app=source->ra_app; 
  burst->dec_app=source->dec_app;
  user->recentre=0; // do not recentre visibilties to burst beam coordinates
  /* default parameteters of the raw data files */
  rfile->have_idx=1;// file header includes the index of the file
  rfile->nfiles=16; // has to be 16 here so that slice_interval etc. is correct
  strcpy(rfile->path,".");
  for(i=0;i<rfile->nfiles;i++)
    sprintf(rfile->fname[i],"%d.dat",i);// NEED TO UPDATE TO CORRECT CONVENTION
  rfile->rec_per_slice=50; // 50 lta records per slice
  // time duration (sec) of each slice
  rfile->t_slice=rfile->rec_per_slice*daspar->lta*user->statime;
  // time interval (sec) between two successive slices
  rfile->slice_interval=rfile->nfiles*rfile->t_slice;
  rfile->t_start[0]=0;
  for(i=1;i<rfile->nfiles;i++) // start time for data in file
    rfile->t_start[i]=rfile->t_start[i-1]+rfile->t_slice;
  rfile->mjd_ref=0.0;// computed when the first slice is read
  user->n_lta=-1;// process all slices (only used when all_data=1)
  user->n_dut1=0.0;// no dut1 correction available by default
#ifdef DO_DUT_CORR
  init_dut1tab(user,bulletinA);
#endif
  //override defaults with actual parameter settings in the binary header
  //file, if specified
  if(strlen(bhdrfile)){
    FILE *fp;
    if((fp=fopen(bhdrfile,"rb"))==NULL)
      {fprintf(stderr,"Unable to open %s\n",bhdrfile); return -1;}
    if((fread(user->corr,sizeof(CorrType),1,fp))!=1)
      {fprintf(stderr,"Error reading corr from %s\n",bhdrfile); return -1;}
    if((fread(user->srec->scan,sizeof(ScanInfoType),1,fp))!=1)
      {fprintf(stderr,"Error reading scan from %s\n",bhdrfile); return -1;}
    // WORK In PROGRESS!! The appropriate parameters are not filled in the
    // binary header. Better to wait till that dust settles before getting
    // into this. JNC 14/June/25
  }
  // over ride earlier settings by parameters given by the user
  if(svuserInp(uparfile,user)) return -1;
  if(user->all_chan) user->channels=corr->daspar.channels;
  if(user->all_data){
    user->channels=corr->daspar.channels/user->nchav;
    source->ch_width *=user->nchav;
  }
  //update burst parameters
  if(user->update_burst)update_burst(user);
  // setup the visibility meta data
  if((user->baselines=init_vispar(user))<=0)
    {fprintf(stderr,"No baselines selected!\n"); return -1;}
  //open the post-correlation beam file
  if(user->postcorr){
    if((user->pcfp=fopen("svfits_postcorr.dat","wb"))==NULL)
      { fprintf(stderr,"Cannot open svfits_postcorr.dat\n"); return -1;}
  }
  if(user->num_threads<0)
    {fprintf(stderr,"Illegal NumThreads %d\n",user->num_threads); return -1;}
#ifdef USE_NOVAS // do coordinate transformation using NOVAS instead of SLA
  init_mat(user,0.0); // initialize matrices for J2000 conversion etc.
#endif
  if(user->epoch<0 && user->recentre)
    {fprintf(stderr,"recentre works only for J2000 epoch"); return -1;}
  return 0;
}
/*
  The start of each visibility slice contains a time_val struct which needs
  to be converted to IST, this is one of the functions used for that. Although
  it says unix_time_to_ist, the f_sec carries the fractional seconds since
  midnight of 01/01/1970, so the returned time is accurate to microseconds
  (i.e. the accuracy of timeval).

  jnc
*/
double unix_time_to_ist(time_t unix_time, double f_sec){
  int    hh,mm,ss;
  double t_ist;
  char   date[256];
  ctime_r(&unix_time,date);//in IST now
  date[strcspn(date,"\n")]='\0';//remove trailing newline,if any
  sscanf(date,"%*s %*s %*s %d:%d:%d %*s",&hh,&mm,&ss);
  t_ist=hh*3600+mm*60+ss+f_sec; //IST seconds
  return t_ist;
}
/*
  computes the reference mjd for the visibility files. The reference mjd is
  the mjd at the earliest midnight IST prior to the timestamp on the first
  slice of the first raw visibility files (which is stored in
  recfile->mjd_ref).

  While computing uvw coordinates etc. the slice time is converted to IST,
  i.e. the time interval since midnight, so converting this to a fractional
  day and adding it to mjd_ref would give the mjd of the visibilty.
  
  jnc may 2025

  fixed bugs in the computation of jd as well as the computation of mjd
  at midnight IST (i.e. is basically UTC, will need correction to get to
  TT)
  jnc 14june 2025
  
*/
double unix_time_to_mjd_ref(time_t unix_time, double f_sec){
  double jd;
  double mjd_ref;

  jd=(unix_time+f_sec)/86400.0+2440587.5; //add JD on 01/01/1970
  mjd_ref=floor(jd-2400000.5);// mjd at midnight UTC
  mjd_ref=mjd_ref-(5.5*3600.0)/86400.0; // mjd at midnight IST
  
  return mjd_ref;
}
/*
  The user structure contains the path to all the raw visibility files
  (assumed to be the same for all files) as well as an array of names of
  the individual files. This function appends the pathname to the raw
  visibility file name of the specified index in the filename array
   and opens it. The calling program should ensure that it is closed when
   done.

   jnc apr 2025
*/
int open_file(SvSelectionType *user,int idx,char *fname,FILE **fp){
  RecFileParType *rfile=&user->recfile;

  strcpy(fname,"");
  if(rfile->path !=NULL)
    {strcpy(fname,rfile->path);fname[strcspn(fname,"\n")]=0;strcat(fname,"/");}
  strcat(fname,rfile->fname[idx]);
  fname[strcspn(fname,"\n")]='\0';//remove trailing newline,if any
  if((*fp=fopen(fname,"rb"))==NULL){
    fprintf(stderr,"Unable to open %s\n",fname); return -1;}
  return 0;
}
/*
  function to read the timestamp embedded in the raw visibility file and
  convert it into IST seconds. The function is expected to be called first
  to read the time stamp of the first slice in the file, which is then stored.
  Subsequent calls to read the timestamp of other slices returns an error if
  the read timestamp differs significantly from the expected one.

  Inputs are the user structure, idx is the index of the file (in the filename
  array) to be read, and slice is the slice number for which the time is
  needed.
  jnc apr 2025

  added computation of mjd_ref for the visibility file,see
  unix_time_to_mjd_ref() for more details.
  jnc may 2025

  idx is now embedded in the file. Modified the code to read it (if
  user->recfile.have_idx is set) and ignore the value supplied to the function

  26/Sep/25 added setting of daspar->mjdref
  
*/
double get_slice_time(SvSelectionType *user, int idx, int slice){
  RecFileParType *rfile=&user->recfile;
  double          slice_interval=rfile->slice_interval;
  int             rec_per_slice=rfile->rec_per_slice;
  double          t_slice=rfile->t_slice;
  unsigned long   timesize=sizeof(struct timeval);
  unsigned long   hdrsize;
  int             have_idx=user->recfile.have_idx;
  unsigned long   recl,off;
  double          start_time;  
  time_t          unix_time;
  double          tiny;
  FILE           *fp;
  char            fname[LINELEN+PATHLEN+1];
  DasParType     *daspar=&user->corr->daspar;
  struct timeval  tv;
  
  if(user->fake_data){// use start time set in init_user() [debug use only]
    start_time=rfile->t_start[idx]+slice*slice_interval;
    return start_time;
  }

  //new format includes the file index in the header
  if(have_idx)hdrsize=timesize+sizeof(int);
  else hdrsize=timesize;
  
  // get the start time of the file
  if(open_file(user,idx,fname,&fp)) return -1;
  recl=user->corr->daspar.baselines*user->corr->daspar.channels*sizeof(float);
  off=slice*(hdrsize+recl*rec_per_slice);
  fseek(fp,off,SEEK_SET);
  if(fread(&tv,timesize,1,fp) != 1)
    {fprintf(stderr,"TimeStmp ReadErr for %s!\n", rfile->fname[idx]);return -1;}
  // read the index to use for timestamp from the file itself. Ignore value
  // supplied while calling the program
  if(have_idx){
    int idx0=idx;
    if(fread(&idx,sizeof(int),1,fp) != 1)
      {fprintf(stderr,"Idx ReadErr for %s!\n", rfile->fname[idx0]);return -1;}
    idx=idx-1;
    if(idx < 0 || idx >= rfile->nfiles){
      fprintf(stderr,"Idx %d RangeErr for %s!\n",idx,rfile->fname[idx0]);
      return -1;
    }
  }
  unix_time=tv.tv_sec; // seconds;
  start_time=unix_time_to_ist(unix_time,tv.tv_usec/1.0e6);// IST
  start_time+=idx*t_slice;//timestamp marks start for for file[0]!
  if(slice==0){
    rfile->t_start[idx]=start_time;
    if(idx==0){//set the mjd_ref for the visibility files
      rfile->mjd_ref=unix_time_to_mjd_ref(unix_time,tv.tv_usec/1.0e6);
      fprintf(user->lfp,"MJD_REF for visibility files = %20.10f\n",
	      rfile->mjd_ref);
      daspar->mjd_ref=rfile->mjd_ref; // USE THIS MJD REF INSTEAD OF USER SUPPLIED
    }
  }else{
    if(start_time<rfile->t_start[idx]){// data crosses midnight
      start_time=start_time+86400.0;// add 1 day, since time is from mjd_ref
      fprintf(user->lfp,"Midnight Crossed idx %d slice %d\n",idx,slice);
    }
    else{//check that timestamp matches to a small fraction of integration time
      tiny=(rfile->t_slice/rfile->rec_per_slice)*1.0e-3;
      if(fabs(start_time-(rfile->t_start[idx]+slice*slice_interval))> tiny){
	fprintf(stderr,"File %d Slice %d TimeStamp Error!\n",idx,slice);
	fprintf(stderr,"Expected %14.6f (sec) Got %14.6f (sec)\n",
		rfile->t_start[idx]+slice*slice_interval,start_time);
	exit(1);
      }
    }
  }
  fclose(fp);

  return start_time;
}
/*
  function that returns the start record as well as the number of contiguous
  records in which the burst signal is present, for the specified raw
  visibility file.

  idx is the index number of the raw visibility file in the array of filenames,
  t_start (contained in user->recfile.t_start[idx], which should be filled
  by an earlier call to get_slice_time()) is the start time of the first slice
  in the file. On return start_rec is the record number of the first record
  in the file which contains the burst, and n_rec is the number of subsequent
  consecutive records which contain the burst. nrec==0 means that this file
  does not contain any burst data. Both start_rec and n_rec are returned in
  the corresponding arrays in user->recfile.

  jnc apr 2025
*/
int get_rec_num(SvSelectionType *user, int idx){
  RecFileParType *rfile=&user->recfile;
  double          t_start=rfile->t_start[idx];
  int            *start_rec=rfile->start_rec+idx;
  int            *n_rec=rfile->n_rec+idx;
  BurstParType   *burst=&user->burst;
  CorrType       *corr=user->corr;
  ScanInfoType   *scan=user->srec->scan;
  double          fb=burst->f/1.0e9;//GHz
  double          burst_width=burst->width;
  double          slice_interval=rfile->slice_interval;
  double          t_slice=rfile->t_slice;
  int             rec_per_slice=rfile->rec_per_slice;
  double          integ=corr->daspar.lta*user->statime;
  double          freq=scan->source.freq[0]/1.0e9; //GHz (first chan)
  double          cw=scan->source.ch_width/1.0e9;//GHz (positive)
  double          t_burst,t_rec,fh,fl;
  int             start_slice,rec_num,r;

  if(scan->source.net_sign[0]>0){fh=freq+corr->daspar.channels*cw;fl=freq;}
  else{fh=freq;fl=freq-corr->daspar.channels*cw;}
  
  // time that the trailing end of the burst arrives at the highest freq
  t_burst=burst->t-burst->int_wd/2.0+K0*burst->DM*(1.0/(fh*fh)-1.0/(fb*fb));
  
  //first slice with data
  start_slice=(int)floor((t_burst-t_start)/slice_interval);
  if(start_slice<0) start_slice=0;
  if(t_burst-(t_start+start_slice*slice_interval) > t_slice){
     //burst starts in the previous slice in some other file
    start_slice++; rec_num=0;
  }else{
    //first record (in selected slice) with data
  rec_num=(int)floor((t_burst-(t_start+start_slice*slice_interval))/integ);
  if(rec_num<0) rec_num=0;
  }
  *start_rec=start_slice*rec_per_slice+rec_num;
  t_rec=t_start+start_slice*slice_interval+rec_num*integ;
  rfile->b_start[idx]=t_rec;//debug use only
  //burst_width is the observed width of the burst, so we do not make
  // any corrrections for intrinsic width, propogation across the band etc.
  for(r=0;t_rec<= t_burst+burst_width;r++){
    if(rec_num == rec_per_slice-1)//start of new slice
      { t_rec +=slice_interval-t_slice; rec_num=0;}
    else
      {t_rec +=integ;rec_num++;}
  }
  *n_rec=r; //n_rec==0 ==> no burst data in this file

  return 0;
}
/*
  function to sort the files in the order of the burst signal arrival time
  in that file. The idea was to avoid having to sort the UV data later 
  into 'TB' order. However this really works only if the burst appears in only
  one slice of a given file. Nonetheless, this function is called to determine
  the section of each file that has to be copied.

  jnc may 2025
*/
int get_file_order(SvSelectionType *user, int *order){
  RecFileParType *rfile=&user->recfile;
  BurstParType   *burst=&user->burst;
  CorrType       *corr=user->corr;
  double          slice_interval=rfile->slice_interval;
  double          t_slice=rfile->t_slice;
  int             rec_per_slice=rfile->rec_per_slice;
  double          integ=corr->daspar.lta*user->statime;
  int             idx,i,j,k,l,index[MaxRecFiles],nfiles;
  double          b_start[MaxRecFiles],min,tmp;

  if(rfile->nfiles>MaxRecFiles){
    fprintf(stderr,"exceeded max allowed files of %d in get_file_order\n",
	    MaxRecFiles);
    return -1;
  }
  // get start_time, start_rec etc for all files
  for(i=0,idx=0;idx<rfile->nfiles;idx++){
    if(get_slice_time(user,idx,0)<0.0)
      {fprintf(stderr,"Could not get time for File %d\n,",idx); return -1;}
    if(get_rec_num(user,idx))
      {fprintf(stderr,"Could not recnum for File %d\n,",idx); return -1;}
    if(rfile->n_rec[idx]>0){
      b_start[i]=rfile->b_start[idx];
      index[i]=idx;
      i++;
    }
  }
  nfiles=i;
  if(nfiles==0){
    fprintf(user->lfp,"No files with burst data!\n");
    return 0; // no files with data
  }

  // order the files in time order using brute force sort
  for(i=0;i<nfiles;i++){
    min=b_start[i];k=i;
    for(j=i+1;j<nfiles;j++){
      if(b_start[j]<=min)
	{min=b_start[j];k=j;}
    }
    order[i]=index[k];
    tmp=b_start[i];
    l=index[i];
    b_start[i]=b_start[k];
    index[i]=index[k];
    b_start[k]=tmp;
    index[k]=l;
  }

  for(i=0;i<nfiles;i++)
    fprintf(user->lfp,"%d %f\n",order[i],rfile->b_start[order[i]]);

  return nfiles;
}
/*
  computes the maximum number of channels over which the signal is spread,
  this is useful for doing buffer memory allocations.
  
  The number of channels decreases with decreasing frequency, i.e. decreases
  with increasing arrival time. However, the first lta record may also not
  cover the maximum number of channels, since it could those parts
  of the pulse which arrived before the recording started. We do an
  measurement over a small (but adhoc chosen) interval hear the start of the
  burst to determine the maximum channel spread.

  This function is not used (in apr 2025, version 0.92), since we have shifted
  to a different logic where we process a complete slice (i.e. 50 contiguous
  records) at a time.

   jnc 15/mar/25

   changed the input parameters to just user, since statime is now inside
   user and not corrpar.
*/
int get_nchan(SvSelectionType *user){
  BurstParType *burst=&user->burst;
  CorrType     *corr=user->corr;
  ScanInfoType *scan=user->srec->scan;
  double        DM=burst->DM;
  double        fb=burst->f/1.0e9; //GHz
  double        tb=0; // bursttime=0, ok since we are just counting channels
  double        wd=burst->int_wd;
  double        freq=scan->source.freq[0]/1.0e9;//GHz (first chan)
  double        cw=scan->source.ch_width/1.0e9; //GHz (positive)
  double        integ=corr->daspar.lta*user->statime;
  int           cs,ce,nc,nrec;
  double        tr1,tr2,fs,fe,fn,fh;
  
  if(scan->source.net_sign[0]>0)fh=freq+corr->daspar.channels*cw;
  else fh=freq;

  nc=0;
  nrec=10;// somewhat adhoc!
  for(tr1=0;tr1<nrec*integ;tr1+=integ){
    tr2=tr1+integ;
    fn=1.0/(fb*fb)+(tr1-tb-wd/2)/(K0*DM);
    fs= (fn>0.0)? 1.0/sqrt(fn):-1.0; // highest frequency containing pulse
    fn=1.0/(fb*fb)+(tr2-tb+wd/2)/(K0*DM);
    fe= (fn>0.0)? 1.0/sqrt(fn):-1.0; //lowest frequency containing pulse
    if(fs>0.0) cs=(int)floor((fh-fs)/cw);
    else cs=-1;
    if(fe>=0.0)ce=(int)floor((fh-fe)/cw);
    else ce=-1;
    if(ce>0 && cs>0)
      if(ce-cs+1>nc) nc=ce-cs+1;
  }
  return nc;
}
/*
  computes the channel range over which the signal is present in a given
  lta record. Only the visibilities corresponding to these channels is
  copied for conversion into FITS.

  jnc 15/mar/25

  changed the input parameters to just user, since statime is now inside
  user and not corrpar.

  jnc june 2025
*/
int get_chan_num(double trec,SvSelectionType *user,int *cs,int *ce){
  BurstParType *burst=&user->burst;
  CorrType     *corr=user->corr;
  ScanInfoType *scan=user->srec->scan;
  double        DM=burst->DM;
  double        tb=burst->t;
  double        fb=burst->f/1.0e9; //GHz
  double        wd=burst->int_wd;
  double        freq=scan->source.freq[0]/1.0e9;//GHz (first chan)
  double        cw=scan->source.ch_width/1.0e9; //GHz (positive)
  double        integ=corr->daspar.lta*user->statime;
  int           c0,c1,c;
  double        tr1,tr2,fs,fe,fn,fh,df;
  
  if(scan->source.net_sign[0]>0)fh=freq+corr->daspar.channels*cw;
  else fh=freq;

  df=cw*scan->source.net_sign[0]; // signed
  tr1=trec; 
  tr2=tr1+integ;
  fn=1.0/(fb*fb)+(tr1-tb-wd/2)/(K0*DM);
  fs= (fn>0.0)? 1.0/sqrt(fn):-1.0; // highest frequency containing pulse
  fn=1.0/(fb*fb)+(tr2-tb+wd/2)/(K0*DM);
  fe= (fn>0.0)? 1.0/sqrt(fn):-1.0; //lowest frequency containing pulse
  if(fs>0.0 && fe> 0.0){
    c0=(int)floor((fs-freq)/df);
    c1=(int)floor((fe-freq)/df);
  }else{
    c0=c1=corr->daspar.channels; // reject all channels
  }
  if(c0>c1){c=c0;c0=c1;c1=c;} // exchange channel limits (df>0)

  // check that the start, end channels and channel range are valid,
  // otherwise reject all channels
  if(c0>0) *cs=c0;else *cs=0;
  if(c1<corr->daspar.channels) *ce=c1;
  else *ce=corr->daspar.channels;
  if(c1<0) // reject all channels
    *cs=*ce=corr->daspar.channels;
  return 0;
}
#define MAT_VEC(A,v,v1){for(i1=0;i1<3;i1++){for(v1[i1]=0,j1=0;j1<3;j1++){v1[i1]+=A[i1][j1]*v[j1];}}}
/*
  rotate the uvw coordinates and correct the visibility phase to the
  burst beam centre in which the burst was found. See init_mat() for the
  initialization of the rotation matrix etc. This effectively changes the
  phase centre from the original one to the burst beam centre.

  It is assumed that the input uvw coordinates are in J2000. The rotation
  is done on the random group format visibilties, which are assumed to be
  in visbuf. The total number of random groups as well as the number of
  channels in each group need to be specified.

  jnc may 2025
*/
void vis_recentre(SvSelectionType *user, char *visbuf, int channels, unsigned long groups){
  UvwParType    *uvwpar;
  ScanInfoType  *scan=user->srec->scan;
  SourceParType *source=&scan->source;
  int            stokes=user->stokes;
  unsigned int   group_size;
  Cmplx3Type    *vis;
  int            g,c,s;
  int            i1,j1; //used in macro MAT_VEC, do *not* reuse
  double         uvw0[3],uvw1[3],re,im,cp,sp,dphs,dphs0;
  double        *lmn=user->lmn_a;
  double         freq0,freq,chwd;

  //frequencies, coordinates etc
  freq0=source->freq[0];
  chwd=source->ch_width*source->net_sign[0];
    
  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes*channels;
#pragma omp parallel for num_threads(user->num_threads) private(g,uvwpar,uvw0,uvw1,vis,c,freq,dphs0,dphs,cp,sp,s,re,im,i1,j1)
  for(g=0;g<groups;g++){
    uvwpar=(UvwParType*)(visbuf+g*group_size);
    uvw0[0]=uvwpar->u; uvw0[1]=uvwpar->v; uvw0[2]=uvwpar->w;
    dphs0=2.0*M_PI*(uvw0[0]*lmn[0]+uvw0[1]*lmn[1]+uvw0[2]*(lmn[2]-1.0));
    MAT_VEC(user->rmat,uvw0,uvw1); //rotate uvw vector to new phase centre
    uvwpar->u=uvw1[0]; uvwpar->v=uvw1[1]; uvwpar->w=uvw1[2];
    vis=(Cmplx3Type*)(visbuf+g*group_size+sizeof(UvwParType));
    // correct the visibility phase for the phase centre shift
    for(c=0;c<channels;c++){
      freq=freq0+c*chwd;//Hz, (uvw are in sec)
      dphs=freq*dphs0;
      cp=cos(dphs);sp=sin(dphs);
      for(s=0;s<stokes;s++){
	if(vis->wt<0) continue;
	re=vis->r*cp+vis->i*sp;
	im=-vis->r*sp+vis->i*cp;
	vis->r=re; vis->i=im;
	vis++;
      }
    }
  }
  return;
}
#undef MAT_VEC
/*
  compute the uvw coordinates. Modified from gvgetUvw in gvfits

  jnc mar25

  moved correction of the precession and nutation to this routine, so
  that the output uvw coordinates from here are in J2000, if user->epoch>0.0
  
  jnc may 2025

  replaced the calculation of the hour angle with the full calculation
  including dut1 and equation of the equinoxes, so that it is now identical
  to what is done in the online software.

  jnc june 2025
*/
int svgetUvw(double tm, SvSelectionType *user, BaseUvwType *uvw)
{ SourceParType *source=&user->srec->scan->source;
  double         mjd_ref=user->recfile.mjd_ref;
  double         lst,dec,ha;
  double         ch,sh, cd,sd ;
  double         bxch,bysh; 
  static double  C = 2.99792458e8 ; /* velocity of light in  m/s */
  int k ;
  
  ha=get_ha(user,tm);  
  dec = source->dec_app;

  ch = cos(ha); sh = sin(ha) ;
  cd = cos(dec) ;  sd = sin(dec) ;

  for (k=0; k<MAX_ANTS; k++)
  { 
    bxch = uvw[k].bx *ch ; /* metres */ 
    bysh = uvw[k].by *sh ;
    uvw[k].u = (uvw[k].bx*sh + uvw[k].by*ch)/C ;          /* sec */
    uvw[k].v = (uvw[k].bz*cd - sd * (bxch - bysh) )/C ;   /* sec */
    uvw[k].w = (cd*(bxch - bysh) + sd*uvw[k].bz)/C  ;     /* sec */
  }

  return 0;
}
/*
  simulate visibilities. Assumes a point source at the location given by
  burst->ra_app,burst->dec_app. Only visibilities for the output baselines
  (i.e. those with antennas in user->antmask) are computed. Visibilities are
  computed at time tm, and for channels between [c0,c1). Output visibilities
  in half float format are returned in rbuf, and the offset appropriate for
  the selected baselines and channels. The calling program should ensure than
  rbuf has enough space allocated.
  
  jnc mar 25

  jnc apr 25 corrected the record structure to not include a timestamp.
             (timestamp is present only once for every recs_per_slice
	     records and is read separately).
*/
int simulate_visibilities(SvSelectionType *user, double tm, char *rbuf,
			  int c0,int c1){
  CorrType        *corr=user->corr;
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  BaseParType     *base=corr->baseline;
  unsigned int     antmask=user->antmask; // antennas in use
  int              channels=corr->daspar.channels;
  int              baselines=corr->daspar.baselines;
  int              stokes=user->stokes;
  int              i,b,c;
  unsigned long    off;
  BaseUvwType      uvw[MAX_ANTS];
  float            freq,u,v,w,l,m,n,re,im,dec,dec0,d_alpha;
  unsigned short  *vis;

  
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}
  if(c0<0 || c0>channels||c1<0 || c1>channels){
      fprintf(stderr,"Illegal ChanRnge %d %d (Legal 0 - %d)\n",c0,c1,
	       channels);return -1;
  }

  //compute uvw coordinates
  for (i=0; i<MAX_ANTS; i++)
  { uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }
  svgetUvw(tm,user,uvw); //units metres/C 

  //compute l,m,n coordinates from AIPS Memo27 Griesen(83,94) SIN projection
  dec=burst->dec_app;
  dec0=source->dec_app;
  d_alpha=burst->ra_app-source->ra_app;
  l=cos(dec)*sin(d_alpha);
  m=sin(dec)*cos(dec0)-cos(dec)*sin(dec0)*cos(d_alpha);
  n=1.0-sqrt(1.0-l*l-m*m);
  //compute visibilities
  for(b=0;b<baselines;b++){
    int a0=base[b].samp[0].ant_id,a1=base[b].samp[1].ant_id;
    int b0=base[b].samp[0].band,b1=base[b].samp[1].band;
    if(a0==a1) continue; // ignore selfs
    if(!(1<<a0 & antmask)) continue; //ignore flagged antennas
    if(!(1<<a1 & antmask)) continue;
    if(b0 !=b1) continue; // ignore cross pol
    for(c=c0;c<c1;c++){
      freq=freq0+c*df;
      u=(uvw[a1].u-uvw[a0].u)*freq;
      v=(uvw[a1].v-uvw[a0].v)*freq;
      w=(uvw[a1].w-uvw[a0].w)*freq;
      re=cos(2.0*M_PI*(u*l+v*m+w*n));
      im=sin(2.0*M_PI*(u*l+v*m+w*n));
      // offset to this channel
      off=(b*channels+c)*sizeof(float); 
      vis=(unsigned short*)(rbuf+off);
      vis[0]=float_to_half(re);
      vis[1]=float_to_half(im);
    }
  }

  return 0;
}
/*
  do a mad based clipping on the amplitudes of the visibilities. The statistics
  are computed on a small radomly selected subset of the visibilities selected
  for conversion to FITS. The parameters idx and slice are for logging purpose
  only to identify what fraction of data from the specfied slice got flagged.
  
  The parameters user->thresh sets the threshold mad above which the visibility
  is flagged. This function also keeps track of the total number of visibilities
  flagged, as well as the number flagged per baseline, both of which are
  reported to the log file.

  jnc apr 2025
*/
int clip(char *visbuf, SvSelectionType *user, int idx, int slice, int groups){
  double         thresh=user->thresh;
  int            stokes=user->stokes;
  Cmplx3Type    *vis;
  UvwParType    *uvwpar;
  long           off=0;
  float         *amp,*med_a,*mad_a,med,mad,a,dum;
  int            sets,set_size,sel_sets,group_size;
  int            i,j,k,flagged=0;
  int            flag_stats[MAX_ANTS][MAX_ANTS];
  int            ant0,ant1,a0,a1;
  
  sel_sets=16;
  set_size=128;
  if((amp=(float*)malloc(set_size*sizeof(float)))==NULL)
    {fprintf(stderr,"Malloc Error\n"); return -1;}
  if((med_a=(float*)malloc(sel_sets*sizeof(float)))==NULL)
    {fprintf(stderr,"Malloc Error\n"); return -1;}
  if((mad_a=(float*)malloc(sel_sets*sizeof(float)))==NULL)
    {fprintf(stderr,"Malloc Error\n"); return -1;}

  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes;


  sets=(int)floor(groups/(1.0*set_size));
  if(sets==0){//only one set
    off=0;
    for(i=0;i<groups;i++){
      off=i*group_size+(i%2)*sizeof(Cmplx3Type); //alternate between RR&LLw
      vis=(Cmplx3Type*)(visbuf+off+sizeof(UvwParType)); //data for this group
      amp[i]=sqrt(vis->r*vis->r+vis->i*vis->i);
    }
    robust_stats(groups,amp,&med,&mad);
  }else{// compute median of medians
    for(i=0;i<sel_sets;i++){
      //randomly pick a set -- we may end up with duplicates here
      j=(random()/RAND_MAX)*sel_sets;
      for(k=0;k<set_size;k++){
	off=(j*set_size+k)*group_size+(k%2)*sizeof(Cmplx3Type);
	vis=(Cmplx3Type*)(visbuf+off+sizeof(UvwParType)); //data for this group
	amp[k]=sqrt(vis->r*vis->r+vis->i*vis->i);
      }
      robust_stats(set_size,amp,&med_a[i],&mad_a[i]);
    }
    robust_stats(sel_sets,med_a,&med,&dum);//median of medians
    robust_stats(sel_sets,mad_a,&mad,&dum);//median of mads
  }

  for(a0=0;a0<MAX_ANTS;a0++)
    for(a1=0;a1<MAX_ANTS;a1++) flag_stats[a0][a1]=0;
  //flag the data
  for(i=0;i<groups;i++){
    off=i*group_size;
    uvwpar=(UvwParType*)(visbuf+off);
    ant0=(int)floor((uvwpar->baseline-257)/256);
    ant1=uvwpar->baseline-257-256*ant0;
    off=off+sizeof(UvwParType);
    vis=(Cmplx3Type*)(visbuf+off);
    a=sqrt(vis[0].r*vis[0].r+vis[0].i*vis[0].i);
    if(fabs(a-med)>thresh*mad)
      {vis[0].wt=-1;flagged++;flag_stats[ant0][ant1]++;}
    a=sqrt(vis[1].r*vis[1].r+vis[1].i*vis[1].i);
    if(fabs(a-med)>thresh*mad)
      {vis[1].wt=-1;flagged++;flag_stats[ant1][ant0]++;}
  }

  for(a0=0;a0<MAX_ANTS;a0++)
    for(a1=a0+1;a1<MAX_ANTS;a1++)
      fprintf(user->lfp,"FLAG:File %2d Slice %4d A0 %3d A1 %3d NFLG_R %8d NFLG_L %8d\n",idx,slice,a0,a1,flag_stats[a0][a1],flag_stats[a1][a0]);

  
  
  
  free(amp);free(med_a);free(mad_a);

  return flagged;
}
/*
  function to compute the "bandpass" and mean "off source" visibility. In
  detail, it returns in the user->bpass structure:(1) the average (over all
  records for which there is no burst signal in the given channel) the real
  and imaginary part separately for each channel of each baseline selected
  for output. This average (arithmentic mean) can be subtracted from the
  signal visibility - the idea being that this would get rid of some RFI
  at least. (2) the average "bandpass" amplitude, which is simply the
  average amplitude of all records in which the signal is not present for
  this channel. It is ensured that the average amplitude bandpass has
  some sensible value in it, even if some channels have no data. The amplitude
  bandpass is normalized using the median value across all channels and
  (3) the total number of visibilities in this slice which contain the signal.
  This is useful for memorly allocation in the calling program.

  It is assumed that rbuf contails all of the data for the given file-index
  and slice.

  jnc apr 2025

  added the option to allow the user to do just one of amplitude bandpass
  (user->do_band True) or off source subtraction (user->do_base True) or
  both (user->do_band True and user->do_base True). Also fixed a bug where
  the bandpass normalization was not happening properly.

  jnc 18apr25
*/
  
int make_bpass(SvSelectionType *user, BpassType *bpass, char *rbuf, int idx, int slice){
  RecFileParType  *rfile=&user->recfile;
  int              rec_per_slice=rfile->rec_per_slice;
  CorrType        *corr=user->corr;
  VisParType      *vispar=&user->vispar;
  VisInfoType     *visinfo=vispar->visinfo;
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  double           integ=corr->daspar.lta*user->statime;
  double           start_time=rfile->t_start[idx]+slice*rfile->slice_interval;
  unsigned long    off,recl;
  int              r,r0,r1,n,b,c,c0,c1,n_vis;
  unsigned short  *in;
  Complex         *off_src;
  float            re,im,median,mad,*abp;
  double           tm;

  n_vis=0; bpass->file_idx=idx; bpass->slice=slice;
  recl=corr->daspar.channels*corr->daspar.baselines*sizeof(float);

  // removed the bypass of this for user->all_chan jnc 5/may/25
  for(r=0;r<rec_per_slice;r++){
    tm=start_time+r*integ;
    get_chan_num(tm,user,&c0,&c1);
    bpass->start_chan[r]=c0;
    bpass->end_chan[r]=c1;
    // compute the number of random groups needed to store the data. Each
    // group has 2 baselines (i.e. stokes)
    if(c1>c0) n_vis +=(baselines/2)*(c1-c0+1);
  }

  if(user->fake_data || !(user->do_band||user->do_base)){//fill nominal values
    for(b=0;b<baselines;b++){
      off_src=bpass->off_src[b]; abp=bpass->abp[b];
      for(c=0;c<channels;c++){
	abp[c]=1.0;
	off_src[c].r=off_src[c].i=0.0;
      }
    }
    return n_vis;
  }

  // for each channel identify the "off source" region, compute the off source
  // signal and the amplitude "bandpass" (i.e. mean over visibilities). In
  // case the user has chosen only one of do_band and do_base the appropriate
  // array is reset to nominal values further below.
#pragma omp parallel for num_threads(user->num_threads) private(c,r0,r1,r,b,off_src,abp,n,off,in,re,im)
  for(c=0;c<channels;c++){
    r0=rec_per_slice;r1=0;
    for(r=0;r<rec_per_slice;r++){
      if(bpass->start_chan[r]<c && bpass->end_chan[r]>c )
	{if(r<r0) r0=r;if(r>r1) r1=r;}
    }
    for(b=0;b<baselines;b++){
      off_src=bpass->off_src[b]; abp=bpass->abp[b];
      off_src[c].r=off_src[c].i=abp[c]=0.0;
      for(n=0,r=0;r<rec_per_slice;r++){
	if(r>=r0 && r<=r1) continue; // burst region
	off=r*recl+visinfo[b].off;
	in=(unsigned short*)(rbuf+off)+2*c;
	re=half_to_float(in[0]);
	im=half_to_float(in[1]);
	if(isfinite(re) && isfinite(im)){
	  off_src[c].r+=re; off_src[c].i+=im;
	  abp[c]+=sqrt(re*re+im*im); n++;
	}
      }
      if(n>0){
	off_src[c].r=off_src[c].r/n;off_src[c].i=off_src[c].i/n;
	abp[c]=abp[c]/n;
      }else{
	off_src[c].r=off_src[c].i=0.0;
	abp[c]=-1.0;// will interpolate over later
	if(b==0)
	  fprintf(user->lfp,"file %d slice %d chan %c No Band/Base\n",
		  idx,slice,c);
      }
    }
  }

  
  if(!user->do_band){//nominal values for bandpass
    for(b=0;b<baselines;b++){
      abp=bpass->abp[b];
      for(c=0;c<channels;c++)
	abp[c]=1.0;
    }
  return n_vis;
  }

  // interpolate bandpass over flagged channels and normalize
  for(b=0;b<baselines;b++){
    abp=bpass->abp[b];
    if(abp[0]<0.0){
      for(c=1;c<channels;c++)
	if(abp[c]>0.0)
	  {abp[0]=abp[c];break;}
    }
    for(c=1;c<channels;c++)
      if(abp[c]<0.0)abp[c]=abp[c-1];
    if(abp[0]<0.0)
      {fprintf(stderr,"All channels flagged, cannot normalize!\n"); return -1;}
    robust_stats(channels,abp,&median,&mad);
    for(c=0;c<channels;c++)
      abp[c]=abp[c]/median;
  }

  if(!user->do_base){// replace off source with nominal values
    for(b=0;b<baselines;b++){
      off_src=bpass->off_src[b];
      for(c=0;c<channels;c++)
	{off_src[c].r=off_src[c].i=0.0;}
    }
    return n_vis;
  }

  bpass->file_idx=idx;
  bpass->slice=slice;
  return n_vis;
}
/*
  read the data for the specified file and slice into rbuf. It is assumed that
  the calling program has allocated sufficient memory for rbuf.

  jnc apr 2025
*/
int read_slice(SvSelectionType *user,int idx, int slice, char *rbuf){
  unsigned long recl;
  unsigned long rec_per_slice=user->recfile.rec_per_slice;
  CorrType     *corr=user->corr;
  unsigned long off;
  unsigned long timesize=sizeof(struct timeval);
  char            fname[LINELEN+PATHLEN+1];
  FILE         *fp;
  // confirm that the slice time is as expected
  if(get_slice_time(user,idx,slice)<0.0)
    return -1;
  recl=corr->daspar.channels*corr->daspar.baselines*sizeof(float);
  off=(unsigned long)slice*(timesize+rec_per_slice*recl)+timesize;//start of data
  if(open_file(user,idx,fname,&fp))return -1;
  if(fseek(fp,off,SEEK_SET)<0)
    {fprintf(stderr,"Error reading %s Slice %d\n",fname,slice); return -1;}
  if(fread(rbuf,recl,rec_per_slice,fp)!=rec_per_slice)
    {fprintf(stderr,"Error reading %s Slice %d\n",fname,slice); return -1;}

  if(0){  // fill with fake data for debug JNC June 2025
    int             r,b,c;
    unsigned short *in;
    int             baselines=corr->daspar.baselines;
    int             channels=corr->daspar.channels;
    for(r=0;r<rec_per_slice;r++){
      for(b=0;b<baselines;b++){
	for(c=0;c<channels;c++){
	  in=(unsigned short*)(rbuf+r*recl+(b*channels+c)*sizeof(float));
	  in[0]=half_to_float(c);in[1]=half_to_float(0);
	}
      }
    }
  }
  return 0;
}

/*
  make a one frequency channel random group, using the channels identified as
  containing the burst. The u,v,w parameters of the group is adjusted so that
  it will be interpreted correctly in the imaging program.

  jnc apr 2025
 */
int make_onechan_group(SvSelectionType *user, BpassType *bpass, int idx, int slice, char *rbuf,
		       int r, double tm, BaseUvwType *uvw,char *obuf,
		       int n_group, unsigned long group_size,int copied0){
  CorrType        *corr=user->corr;
  VisParType      *vispar=&user->vispar;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  double           epoch,epoch1=user->epoch;
  int              c,c0,c1,copied,b,b1;
  unsigned long    recl,off;
  char            *rbuf1;
  double           JD,date2;
  int              date1;
  int              baselines=user->baselines;//output baselines
  int              channels=corr->daspar.channels;//input channels
  double           integ=corr->daspar.lta*user->statime;
  Cmplx3Type      *out;
  unsigned short  *in;
  UvwParType      *uvwpar;
  Complex         *off_src;
  float           *abp,freq;
  SvSelectionType  user1; //local copy

  copied=copied0;//groups already copied
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  c0=bpass->start_chan[r];c1=bpass->end_chan[r];
  if(user->fake_data){
    //simulate the visibilities
    if(simulate_visibilities(user,tm,rbuf,c0,c1)) return -1;
    rbuf1=rbuf;
  }else{ // go to this record
    rbuf1=rbuf+r*recl;
  }
  if(user->do_log>10)
    fprintf(user->lfp,"COPY:File %d Slice %d Rec %d Time %12.6f Chan %d - %d\n",
	    idx,slice,r,tm,c0,c1);
  if(c0<0 || c0>=channels) return 0; // record does not contain burst
  svgetUvw(tm,user,uvw);
#ifdef USE_NOVAS
  // initialize the rotation matrix to J2000 etc. Copy over to user1 to be
  // thread safe.
  memcpy(&user1,user,sizeof(SvSelectionType));
  init_mat(&user1,tm); 
#endif
  JD = user->recfile.mjd_ref + 2400000.5 + tm/86400 ;
  date1 = (int) JD ;
  date2 = JD - date1+user->iatutc/86400 ; 
  for(b=0;b<baselines;b+=vispar->vis_sets){//vis_sets==2
    VisInfoType *vinfo=vispar->visinfo;
    int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
    for(c=c0;c<c1;c++){ // limits for this timestamp
      if(c>=channels) return copied-copied0;
      freq=freq0+c*df; // frequency of this channel
      if(copied ==n_group){
	fprintf(stderr,"Exceeded the estimated number of groups [%d]!\n",
		n_group);
	return -1;
      }
      off=copied*group_size;
      uvwpar=(UvwParType*)(obuf+off);//uvw for this group
      out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType));//data for this group
      // Note off_src, abp are set to (0,1) in case the user does not
      // want to do calibration
      for(b1=b;b1<b+vispar->vis_sets;b1++){
	off_src=bpass->off_src[b1];
	abp=bpass->abp[b1];
	in=(unsigned short*)(rbuf1+vinfo[b1].off)+2*c;//input visibility
	out->r=half_to_float(in[0]); out->i=half_to_float(in[1]);
	if(isfinite(out->r) && isfinite(out->i)){
	  out->r=(out->r-off_src[c].r)/abp[c];
	  out->i=(out->i-off_src[c].i)/abp[c];
	  out->wt=1.0;
	  if(vinfo[b1].flip)out->i=-out->i;
	}else
	  out->wt=-1.0;
	if(vinfo[b1].drop)out->wt=-1.0;
	out++;//now copy same chan for next stokes
      }
      // scale uv co-ordinates so that the values come out right inside of
      // AIPS/CASA. In the imaging programs the number here is multiplied
      // by the frequency of the channel inorder to get the uv coordinate
      // for that channel). In the FITS file, all of the data is put in
      // the first channel, i.e. at frequency freq0. The scaling here
      // corrects for this.
      uvwpar->u = (uvw[ant1].u-uvw[ant0].u)*(freq/freq0);
      uvwpar->v = (uvw[ant1].v-uvw[ant0].v)*(freq/freq0);
      uvwpar->w = (uvw[ant1].w-uvw[ant0].w)*(freq/freq0);
      if(epoch>0.0){//rotate to J2000
#ifdef USE_NOVAS
	novas_prenut_vis(&user1,uvwpar,user->recfile.mjd_ref+tm/86400.0);//user1 to be thread safe;
#else
	sla_prenut_vis(uvwpar,user->recfile.mjd_ref,source->ra_app,
		    source->dec_app,2000.0);
#endif
      }
      uvwpar->date1 = date1 ;
      uvwpar->date2 = date2;
      uvwpar->baseline = ant0*256 + ant1 + 257 ;
      uvwpar->su = 1; // only one source in file
      uvwpar->fq = 1 ; // only one frequency id in file
      copied++;
    }
  }
  return copied-copied0; // number of groups copied in this call
}
/*
  make a postcorrelation beam, this is useful for cross checks. The
  post-correlation beam amplitudes (per channel) are written out, and can
  be plotted (using for example plot_postcorr.py, which also has an option
  for reading in the svfits.log file and over plotting the region identified
  as containing the burst signal).

  jnc apr 2025
*/
int make_postcorr_beam(SvSelectionType *user, BpassType *bpass, char *rbuf,int r, double tm,
		       BaseUvwType *uvw){
  CorrType        *corr=user->corr;
  VisParType      *vispar=&user->vispar;
  VisInfoType     *vinfo=vispar->visinfo;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0],freq;
  double           df=source->ch_width*source->net_sign[0];//signed
  double           cp,sp,re1,im1,dphs[MAX_ANTS];
  double          *lmn_a=user->lmn_a;
  int              k,c,b,p;
  unsigned long    recl,off;
  char            *rbuf1;
  int              baselines=user->baselines;//output baselines
  int              channels=corr->daspar.channels;//input channels
  Cmplx3Type      *out;
  unsigned short  *in;
  Complex         *off_src;
  float           *abp;   //MAX_CHANS in newcorr.h is 8*4096
  float            phs,re,im,re0,im0,bm[MAX_CHANS/8];
  int              ant0,ant1,drop_csq=user->drop_csq;
  static int       wrote_hdr=0; //not thread safe
  unsigned int     antmask=user->antmask;
  //structure is double aligned
  struct postcorr_hdr{double mjd,int_wd,DM,f_start,f_end,integ;long channels;};
  struct postcorr_hdr  phdr;

  if(user->recentre){//compute per antenna phase
    for(k=0;k<MAX_ANTS;k++)
      dphs[k]=2.0*M_PI*(uvw[k].u*lmn_a[0]+uvw[k].v*lmn_a[1]+
			uvw[k].w*(lmn_a[2]-1));
  }else{
    for(k=0;k<MAX_ANTS;k++)dphs[k]=0.0;
  }
   
  // write out some header information on first invocation
  if(!wrote_hdr){
    SourceParType *source=&user->srec->scan->source;
    double         df=source->ch_width*source->net_sign[0];//signed
    double          integ=corr->daspar.lta*user->statime;
    phdr.mjd=user->burst.mjd;
    phdr.int_wd=user->burst.int_wd;
    phdr.DM=user->burst.DM;
    phdr.channels=corr->daspar.channels;
    phdr.f_start=source->freq[0];
    phdr.f_end=source->freq[0]+phdr.channels*df;
    phdr.integ=integ;
    if(fwrite(&phdr,sizeof(struct postcorr_hdr),1,user->pcfp) !=1)
      {perror(NULL);fprintf(stderr,"Error writing hdr to beam output file\n"); return -1;}
    wrote_hdr=1;
  }

  // compute the post-correlation beam for this record
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  rbuf1=rbuf+r*recl;
#pragma omp parallel for num_threads(user->num_threads) private(c,freq,re0,im0,p,b,ant0,ant1,off_src,abp,in,re,im,re1,im1,cp,sp)
  for(c=0;c<channels;c++){
    freq=freq0+c*df;
    re0=im0=bm[c]=0.0;
    for(p=0;p<user->stokes;p++){
      for(b=p;b<baselines;b+=vispar->vis_sets){//vis_sets==user->stokes
	ant0=vinfo[b].ant0; ant1=vinfo[b].ant1;
	if(!(1<<ant0&antmask)||!(1<<ant1&antmask)) continue;
	if(drop_csq && ant0 < 11 && ant1 < 11) continue;
	if(vinfo[b].drop) continue;
	off_src=bpass->off_src[b]; abp=bpass->abp[b];
	in=(unsigned short*)(rbuf1+vinfo[b].off)+2*c;//input visibility
	re=half_to_float(in[0]); im=half_to_float(in[1]);
	if(!(isfinite(re) && isfinite(im))) continue;
	re=(re-off_src[c].r)/abp[c];
	im=(im-off_src[c].i)/abp[c];
	phs=(dphs[ant1]-dphs[ant0])*freq;
	if(vinfo[b].flip){im=-im;phs=-phs;}
	cp=cos(phs);sp=sin(phs);
	re1=re*cp-im*sp;
	im1=re*sp+im*cp;
	re0+=re1;im0+=im1;
      }
      bm[c]+=sqrt(re0*re0+im0*im0);
    }
  }

  //write out the timestamp and post-correlation beam 
  if(fwrite(&tm,sizeof(double),1,user->pcfp) !=1)
    {perror(NULL);fprintf(stderr,"Error writing to beam output file\n"); return -1;}
  if(fwrite(bm,sizeof(float),channels,user->pcfp) !=channels)
    {perror(NULL);fprintf(stderr,"Error writing to beam output file\n"); return -1;}


  return 0;
}
/*
  make a standard multi frequency channel random group using all input channels.
  also makes a post-correlation beam, if requested. This function is called
  whenever one wants to make a post_correlation beam, even if one does not
  want to write a multi-channel output FITS file. In that case the function
  returns after making the post-correlation beam.
  
  jnc apr 2025

  modified to flag the channels that do not contain the burst signal
  jnc may 2025

  fixed the output visibility location, earlier all channels for a given
  stokes were contiguous, which is not as per the FITS format.
  jnc 14june25
 */
int make_allchan_group(SvSelectionType *user, BpassType *bpass, int idx, int slice, char *rbuf,
		       int r, double tm, BaseUvwType *uvw,char *obuf,
		       int n_group, unsigned long group_size,unsigned
		       long copied0){
  CorrType        *corr=user->corr;
  VisParType      *vispar=&user->vispar;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  RecFileParType  *rfile=&user->recfile;
  int              rec_per_slice=rfile->rec_per_slice;
  double           slice_interval=rfile->slice_interval;
  double           integ=corr->daspar.lta*user->statime;
  double           epoch,epoch1=user->epoch;
  int              c,c0,c1,b,b1,s;
  unsigned long    copied;
  unsigned long    recl,off;
  char            *rbuf1;
  double           JD,date2;
  int              date1;
  int              baselines=user->baselines;//output baselines
  int              channels=corr->daspar.channels;//input channels
  Cmplx3Type      *out;
  unsigned short  *in;
  UvwParType      *uvwpar;
  Complex         *off_src;
  float           *abp;
  SvSelectionType user1;
  if(user->postcorr)
    if(make_postcorr_beam(user,bpass,rbuf,r,tm,uvw)) return -1;

  if(!user->all_chan) return 0;//user only wanted a post-correlation beam
  
  copied=copied0;//number of groups already copied
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  c0=bpass->start_chan[r];c1=bpass->end_chan[r];
  rbuf1=rbuf+r*recl;
  // SET TM AS PER CURRENT UNDERSTANDING JNC 2/OCT/25
  tm=rfile->t_start[idx]+slice*rfile->slice_interval+(r+0.5)*integ; //middle of integ
  svgetUvw(tm,user,uvw);
#ifdef USE_NOVAS
  memcpy(&user1,user,sizeof(SvSelectionType)); // just to be thread safe
  init_mat(&user1,tm);// initialize the rotation matrix to J2000 etc.
#endif
  JD = user->recfile.mjd_ref + 2400000.5 + tm/86400.0;
  date1 = (int) JD ;
  date2 = JD - date1+user->iatutc/86400.0 ; 

  for(b=0;b<baselines;b+=vispar->vis_sets){//vis_sets==2
    VisInfoType *vinfo=vispar->visinfo;
    off=copied*group_size;
    uvwpar=(UvwParType*)(obuf+off);//uvw for this group
    if(copied ==n_group){
      fprintf(stderr,"Exceeded the estimated number of groups [%d]!\n",
	      n_group);
      return -1;
    }
    int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
    for(s=0,b1=b;b1<b+vispar->vis_sets;b1++,s++){
      //output location for this group,stokes
      out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType))+s;
      // Note off_src, abp are set to (0,1) in case the user does not
      // want to do calibration
      off_src=bpass->off_src[b1];
      abp=bpass->abp[b1];
      for(c=0;c<channels;c++){
	in=(unsigned short*)(rbuf1+vinfo[b1].off)+2*c;//input visibility
	out->r=half_to_float(in[0]); out->i=half_to_float(in[1]);
	if(isfinite(out->r) && isfinite(out->i)){
	  out->r=(out->r-off_src[c].r)/abp[c];
	  out->i=(out->i-off_src[c].i)/abp[c];
	  out->wt=1.0;
	  if(vinfo[b1].flip)out->i=-out->i;
	}else{out->wt=-1.0;}
	if(vinfo[b1].drop)out->wt=-1.0;
	if(c<c0||c>c1) out->wt=-1.0; // flag regions outside burst
	out+=vispar->vis_sets; // same as stokes
      }
    }
    uvwpar->u = (uvw[ant1].u-uvw[ant0].u);
    uvwpar->v = (uvw[ant1].v-uvw[ant0].v);
    uvwpar->w = (uvw[ant1].w-uvw[ant0].w);
    if(epoch>0.0){
#ifdef USE_NOVAS
      novas_prenut_vis(&user1,uvwpar,user->recfile.mjd_ref+tm/86400.0);
#else
      sla_prenut_vis(uvwpar,user->recfile.mjd_ref,source->ra_app,
		     source->dec_app,2000.0);
#endif
    }
    uvwpar->date1 = date1 ;
    uvwpar->date2 = date2;
    uvwpar->baseline = ant0*256 + ant1 + 257 ;
    uvwpar->su = 1; // only one source in file
    uvwpar->fq = 1 ; // only one frequency id in file
    copied++;
  }
  return copied-copied0; // number of groups copied in this call
}
/*
  function to read in a single SPOTLIGHT raw visibility file (specified by
  idx) and extract only the data containg the burst in set of FITS random
  groups.

  Each channel of each set of baselines (i.e. all baselines corresponding
  to the same antenna pair but with differing stokes) is treated as a
  separate random group, with the UV coordinates in the random parameters
  appropriately scaled so as to give the correct values inside the AIPS/CASA
  imaging package. The UV coordinates are also rotated to J2000 from the
  apparent coordinates (which is the frame in which the visibilities were
  computed).

  The return value of the program is the number of random groups copied, or
  -1 in case of error. The random groups themselves are in outbuf.

  Rawdump visibility files are assumed to have the following format
  TV,b0s0c0,b0s0c1,...b0s1c0,...b1s0c0...
  where TV is a struct timeval giving the timestamp, b is the baseline
  between a pair of antennas, s is the stokes and c is the channel number.
  The real and imaginary parts of the visibility are assumed to be half
  float each. See init_vispar for more details on how the input baseline
  order is mapped onto the output random group order.

  jnc mar 2025

  did a major reorganization of the logic to process one slice at a time.
  (an intermediate version, archived but with the code completely replaced
  here copied a fixed chunk of data at at time, with the chunk being fixed
  by the user and unrelated to the file slice size). Also added the option
  of subtracting the mean and doing amplitude "bandpass" calibration (i.e.
  dividing by the mean amplitude of the off source visibilities), and also
  integrated the code for replacing observed visiblities with simulated ones
  (i.e. with user->fake_data==1) into this function.

  jnc apr 2025
*/
int copy_vis(SvSelectionType *user, BpassType *bpass, int idx, int slice,
	       int start_rec, int n_rec, char *rbuf,char **outbuf){
  UvwParType      *uvwpar;
  RecFileParType  *rfile=&user->recfile;
  int              rec_per_slice=rfile->rec_per_slice;
  double           slice_interval=rfile->slice_interval;
  CorrType        *corr=user->corr;
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  int              group_size,c,n_group,rec0;
  unsigned long    copied;
  int              r,i;
  double           tm;
  double           integ=corr->daspar.lta*user->statime;
  unsigned long    off,bufsize,recl,timesize=sizeof(struct timeval);
  BaseUvwType      uvw[MAX_ANTS];
  char            *rbuf1,*obuf;
  FILE            *lfp=user->lfp,*fp;
  
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}
  if(idx<0 || idx>rfile->nfiles){
    fprintf(stderr,"Illegal raw file number %d [legal %d - %d]\n",
	    idx,0,rfile->nfiles);
    return -1;
  }

  // copy over the X,Y,Z co-ordinates into a local array for computing
  // u,v,w later
  for (i=0; i<MAX_ANTS; i++)
  { uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }

  // compute the off source spectra, the mean amplitude "bandpass" and
  // the number of visibilities (2*groups)to copy. Nominal values
  // (off_src=0, bandpass=1) are returned in case the user does not want
  // to do bandpass calibration.
  n_group=make_bpass(user,bpass,rbuf,idx,slice);
  if(!user->all_chan){
    if(n_group==0){
      fprintf(lfp,"No visibilities copied from %s slice %d\n",
	      rfile->fname[idx],slice);
      return 0;
    }
    // each random group contains data for all stokes and one channel 
    group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes;
    bufsize=n_group*group_size;
    obuf=(*outbuf=(char*)malloc(bufsize));//mem returned to calling prog
    if(obuf==NULL) {fprintf(stderr,"Malloc Error\n");return -1;}
  }else{// allocate memory for copying all channels
    //2 baselines [i.e. polarizations] per group
    n_group=rec_per_slice*baselines/2; 
    group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes*channels;
    bufsize=n_group*group_size;
    obuf=(*outbuf=(char*)malloc(bufsize));//mem returned to calling prog
    if(obuf==NULL) {fprintf(stderr,"Malloc Error\n");return -1;}
  }
    
  // set start record along with its time stamp
  rec0=start_rec-slice*rec_per_slice;
  if(rec0<0) rec0=slice*rec_per_slice;//in the second or later slice with signal
  tm= rfile->t_start[idx]+slice*rfile->slice_interval+rec0*integ;
  tm+=user->timestamp_off; // middle of the record
  // convert selected data into random groups
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  for(copied=0,r=rec0;r<rec_per_slice;r++){
    if(user->all_chan||user->postcorr)
      c=make_allchan_group(user,bpass,idx,slice,rbuf,r,tm,uvw,obuf,n_group,
			 group_size,copied);
    if(!user->all_chan)
      c=make_onechan_group(user,bpass,idx,slice,rbuf,r,tm,uvw,obuf,n_group,
			 group_size,copied);
    if(c<0){
      fprintf(stderr,"Error copying visibilities from %s Slice %d\n",
	      rfile->fname[idx],slice);
      return -1;
    }else copied+=c;
    if (r+slice*rec_per_slice==start_rec+n_rec)//copied all data
      break;
    tm=tm+integ;
  }
  if(user->recentre)// shift phase centre
    vis_recentre(user, obuf,1, copied);
  return copied; // total number of groups
}
enum{SampleSize=1000};
typedef struct vis_sel_type{int r,c;} VisSelType;
int approx_stats(char *rbuf,unsigned int off, int recl,VisSelType *selvis,
		 int nsamp,float *med, float *mad){
  int             v,n;
  float           re,im,vdata[SampleSize];
  unsigned short *in;

  for(n=0,v=0;v<nsamp;v++){
    in=(unsigned short*)(rbuf+selvis[v].r*recl+off)+2*selvis[v].c;
    re=half_to_float(in[0]);im=half_to_float(in[1]);
    if(isfinite(re) && isfinite(im))
      {vdata[n]=sqrt(re*re+im*im);n++;}
  }
  if(robust_stats(n,vdata,med,mad)<0) return -1;
  else return 0;
}
/*
  function to average the visiblties (over all records in the given slice)
  and return them in a random group structure.
  idx is the index of the input raw visibility file, slice is the slice number
  to operate on, rbuf contains all the raw visibility data for this slice, and
  outbuf (memory allocated here and should be freed in the calling program)
  contains the averaged data in random group format. Returns the number of
  random groups created, -1 in case of failure.

  jnc may 2025

  fixed the output visibility location, earlier all channels for a given
  stokes were contiguous, which is not as per the FITS format.
  jnc 14june25

*/
int avg_vis(SvSelectionType *user, int idx, int slice, char *rbuf,
	    char *outbuf){
  UvwParType      *uvwpar;
  RecFileParType  *rfile=&user->recfile;
  int              rec_per_slice=rfile->rec_per_slice;
  double           slice_interval=rfile->slice_interval;
  CorrType        *corr=user->corr;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  VisParType      *vispar=&user->vispar;
  int              channels=user->channels;
  int              nchav=user->nchav;
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  unsigned int     group_size,flagged;
  unsigned long    n_group;
  int              r,i,b,b1,c,c1,n,s;
  double           tm,epoch,epoch1,JD,date2;
  int              date1;
  float            re,im,amp,med,mad;
  double           integ=corr->daspar.lta*user->statime;
  unsigned long    off,bufsize,recl,timesize=sizeof(struct timeval);
  BaseUvwType      uvw[MAX_ANTS];
  char            *obuf=outbuf;
  unsigned short  *in;
  Cmplx3Type      *out;
  FILE            *lfp=user->lfp,*fp;
  VisSelType       selvis[SampleSize];
  SvSelectionType  user1; // local copy
    

    
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}
  if(idx<0 || idx>rfile->nfiles){
    fprintf(stderr,"Illegal raw file number %d [legal %d - %d]\n",
	    idx,0,rfile->nfiles);
    return -1;
  }

  // copy over the X,Y,Z co-ordinates into a local array for computing
  // u,v,w later
  for (i=0; i<MAX_ANTS; i++){
    uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }

  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*(channels*stokes);
  n_group=baselines/stokes;

  //set the timestamp to the middle of the slice 
  tm= rfile->t_start[idx]+slice*rfile->slice_interval;
  tm+=user->timestamp_off+rec_per_slice*(integ/2); // middle of the slice
  svgetUvw(tm,user,uvw);
#ifdef USE_NOVAS
  memcpy(&user1,user,sizeof(SvSelectionType)); // just to be thread safe
  init_mat(&user1,tm);// initialize the rotation matrix to J2000 etc.
#endif
  JD = user->recfile.mjd_ref + 2400000.5 + tm/86400.0;
  date1 = (int) JD ;
  date2 = JD - date1+user->iatutc/86400 ; 

  if(user->do_flag){
    // select some random visibilities to compute statstics for flagging
    srand(time(NULL));
    for(r=0;r<1000;r++){
      selvis[r].r=rand()%rec_per_slice;
      selvis[r].c=rand()%user->corr->daspar.channels;
    }
  }

  // convert selected data into random groups
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  //flagged=0;

  // average over records and channels, the flagged count will not be
  // correct for multi-threaded computation
#pragma omp parallel for num_threads(user->num_threads) private(b,off,uvwpar,out,b1,med,mad,c,n,s,c1,r,in,re,im,amp)
  for(b=0;b<baselines;b+=stokes){//(stokes=2==vispar->vis_sets)
    VisInfoType *vinfo=vispar->visinfo;
    int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
    off=b/stokes*group_size;
    uvwpar=(UvwParType*)(obuf+off);//uvw for this group
    for(s=0,b1=b;b1<b+stokes;b1++,s++){//all baselines in group(i.e. all stokes)
      // start of output data for this stokes in this group
      out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType))+s;
      if(user->do_flag)
	if(approx_stats(rbuf,vinfo[b1].off,recl,selvis,SampleSize,&med,&mad)<0)
	  exit(-1);
      for(c=0;c<corr->daspar.channels;c+=nchav){
	out->r=out->i=0.0;out->wt=-1.0;//default flagged
	for(n=0,c1=c;c1<c+nchav;c1++){//average over channels and records
	  for(r=0;r<rec_per_slice;r++){
	    in=(unsigned short*)(rbuf+r*recl+vinfo[b1].off)+2*c1;//input vis
	    re=half_to_float(in[0]); im=half_to_float(in[1]);
	    if(isfinite(re) && isfinite(im)){
	      if(user->do_flag){
		amp=sqrt(re*re+im*im);
		if(fabs(amp-med)>user->thresh*mad) continue;
	      }
	      out->r+=re; out->i+=im;n++;
	    }
	  }
	}
	//flagged+=(nchav*rec_per_slice-n);
	if(n>0){out->r/=n;out->i/=n;out->wt=(1.0*n)/(nchav*rec_per_slice);}
	if(!vinfo[b1].drop) out->wt=1.0;
	if(vinfo[b1].flip)out->i=-out->i;
	out+=stokes;//RR,LL for a given channel are adjacent
      }
    }
    uvwpar->u = (uvw[ant1].u-uvw[ant0].u);
    uvwpar->v = (uvw[ant1].v-uvw[ant0].v);
    uvwpar->w = (uvw[ant1].w-uvw[ant0].w);
    if(user->epoch>0){ //rotate to J2000
#ifdef USE_NOVAS
      novas_prenut_vis(&user1,uvwpar,user->recfile.mjd_ref+tm/86400.0); //user1 to be thread safe
#else
      sla_prenut_vis(uvwpar,user->recfile.mjd_ref,source->ra_app,
		  source->dec_app,2000.0);
#endif
    }
    uvwpar->date1 = date1 ;
    uvwpar->date2 = date2;
    uvwpar->baseline = ant0*256 + ant1 + 257 ;
    uvwpar->su = 1; // only one source in file
    uvwpar->fq = 1 ; // only one frequency id in file
  }

  if(user->do_log>10){
    fprintf(user->lfp,"AVG_VIS: File %d Slice %d Time %16.7f %16d %26.10lf\n",
	    idx,slice,tm,date1,date2);
    fprintf(user->lfp,"Flagged %d out of %d visibilities\n",flagged,
	    baselines*channels*rec_per_slice);
  }
  if(user->recentre)// move phase centre to burst beam centre
    vis_recentre(&user1, obuf,user->channels, n_group); //user1 to be thread safe
  return n_group; // total number of groups
}
