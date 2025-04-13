#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<time.h>

#include"fitsio.h"
#include "newcorr.h"
#include "svio.h"

enum{Hdrlen=81920,Len=80};
static char Buf[Hdrlen], *Bufp ;// local buffer for reading hdr files
const char cont_char = '*', begin_char = '{', end_char = '}';

// The enum and the char array need to match exactly, and are used to
// parse the input parameter file.
enum {NFiles,Path,InputFile,FitsFile,AntMask,AllChan,IATUTC,Epoch,ObsMjd,
      FreqSet,CoordType,RaApp,DecApp,RaMean,DecMean,StokesType,BurstName,
      BurstMJD,BurstTime,BurstDT,BurstIntWd,BurstWd,BurstDM,BurstDDM,BurstFreq,
      BurstBmId,BurstRa,BurstDec,UpdateBurst,DoFlag,Thresh,DoBand,
      SvVars} SvVarType ;
char *SvVarname[SvVars] =
  {"NFILE","PATH","INPUT","FITS","ANTMASK","ALL_CHAN","IATUTC","EPOCH",
   "OBS_MJD","FREQ_SET","COORD_TYPE","RA_APP","DEC_APP","RA_MEAN","DEC_MEAN",
   "STOKES_TYPE","BURST_NAME","BURST_MJD","BURST_TIME","BURST_DT","BURST_INTWD",
   "BURST_WIDTH","BURST_DM","BURST_DDM","BURST_FREQ","BURST_BM_ID",
   "BURST_RA","BURST_DEC","UPDATE_BURST","DO_FLAG","THRESH","DO_BAND"};
static double K0=4.15e-3; // DM constant in seconds

static float RefFreq,ChanWidth; // NEED TO SET THESE!!!

extern void sla_nut_(double *mjd, double *a); 
extern void sla_prec_(double *epoch0,double *epoch1, double *a); 
extern double sla_epj_(double *mjd);
extern void   sla_amp_(double *ra_app, double *dec_app, double *mjd, double *epoch1,
		       double *ra_mean, double *dec_mean);

#define TINY 1.0e-9
void prenut(UvwParType *uv,double mjd,double ra_app,double dec_app,
	    double epoch1)
{ double  epoch,v[3],v1[3],a[9],nm[3][3],pm[3][3],t,ra_mean,dec_mean;
  static double  p[3][3],p1[3][3],p2[3][3],p3[3][3],rm[3][3];
  int    i,j,k;
  static int new_epoch,new_source;
  static double r=-4*M_PI,d=-4*M_PI,m=0.0,e=0.0;

  new_epoch=new_source=0;
  if(fabs(m-mjd) >TINY ||fabs(e-epoch1) > TINY)
  { new_epoch=1; new_source=1;
    m=mjd;e=epoch1;r=ra_app;d=dec_app;
  }else
  { if(fabs(r-ra_app)>TINY && fabs(d-dec_app)>TINY)
    { new_source=1; r=ra_app;d=dec_app;}
  }
  

  if(new_epoch)
  { sla_nut_(&mjd,a); 
    for(i=0;i<3;i++)   /* nutation to mean on given mjd */
      for(j=0;j<3;j++)nm[j][i]=a[3*i+j];
    epoch=sla_epj_(&mjd); 
    sla_prec_(&epoch,&epoch1,a); 
    for(i=0;i<3;i++) /* precession of mean between epochs */
      for(j=0;j<3;j++)pm[j][i]=a[3*i+j];
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
      { rm[i][j]=0.0;  /* combined rotation matrix */
        for(k=0;k<3;k++)rm[i][j] +=pm[i][k]*nm[j][k];
      }
  }

  if(new_source)
  { /* rotate from source(dec) to equitorial plane */
    t=(M_PI/2-dec_app);
    p[0][0]=1.0; p[0][1]=0.0;     p[0][2]= 0.0;
    p[1][0]=0.0; p[1][1]=cos(t);  p[1][2]=-sin(t);
    p[2][0]=0.0; p[2][1]=sin(t);  p[2][2]= cos(t);
    
    /* rotate from source(ra) to ra=0.0 */
    t=(M_PI/2.0+ra_app);
    p1[0][0]=cos(t); p1[0][1]=-sin(t); p1[0][2]=0.0;
    p1[1][0]=sin(t); p1[1][1]= cos(t); p1[1][2]=0.0;
    p1[2][0]=0.0;    p1[2][1]= 0.0;    p1[2][2]=1.0;

    sla_amp_(&ra_app,&dec_app,&mjd,&epoch1,&ra_mean,&dec_mean);

    /* rotate from source(dec) to equitorial plane */
    t=(M_PI/2-dec_mean);
    p2[0][0]=1.0; p2[0][1]=0.0;     p2[0][2]= 0.0;
    p2[1][0]=0.0; p2[1][1]=cos(t);  p2[1][2]=-sin(t);
    p2[2][0]=0.0; p2[2][1]=sin(t);  p2[2][2]= cos(t);

    /* rotate from source(ra) to ra=0.0 */
    t=(M_PI/2.0+ra_mean);
    p3[0][0]=cos(t); p3[0][1]=-sin(t); p3[0][2]=0.0;
    p3[1][0]=sin(t); p3[1][1]= cos(t); p3[1][2]=0.0;
    p3[2][0]=0.0;    p3[2][1]= 0.0;    p3[2][2]=1.0;
  }

  
  v[0]=uv->u;v[1]=uv->v;v[2]=uv->w;
  /* compute co-ordinates in the equitorial system */
  for(i=0;i<3;i++)
  { v1[i]=0.0;
    for(j=0;j<3;j++)v1[i] += p[i][j]*v[j];
  }
  for(i=0;i<3;i++)
  { v[i]=0.0;
    for(j=0;j<3;j++)v[i] += p1[i][j]*v1[j];
  }
  /* equitorial co-ordinates at the new epoch */
  for(i=0;i<3;i++)
  { v1[i]=0.0;
    for(j=0;j<3;j++)v1[i] += rm[i][j]*v[j];
  }
  /* compute interferometric co-ordinates at the new epoch */
  for(i=0;i<3;i++)
  { v[i]=0.0;
    for(j=0;j<3;j++)v[i] += p3[j][i]*v1[j];
  }
  for(i=0;i<3;i++)
  { v1[i]=0.0;
    for(j=0;j<3;j++)v1[i] += p2[j][i]*v[j];
  }
  uv->u=v1[0];uv->v=v1[1];  uv->w=v1[2];

  return;
}
#undef TINY

float svVersion(void){
  return 1.00;
}
int svuserInp (char *filename, SvSelectionType *user ){
  CorrType      *corr=user->corr;
  DasParType    *daspar=&user->corr->daspar;
  SourceParType *source=&user->srec->scan->source;
  char           str[512],key_str[32] ;
  int            val[512];
  FILE          *fp = fopen(filename, "rt") ;
  int            var, np, k ;
 
  while ( fgets(str, 500, fp) )
  { char *p = str ;
    while (*p && (*p == ' '))p++ ;
    if (*p == 0)continue ;
    if ( (*p == '#') || (*p == '*') )continue ;
    key_str[0] = 0 ;
    sscanf(p, "%s",key_str) ;
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
      case AntMask  :sscanf(p,"%u",&user->antmask);break;
      case AllChan  :sscanf(p,"%d",&user->all_chan);break;
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
      case CoordType:strcpy(user->coord_type,p);break;
      case RaApp    :sscanf(p,"%lf",&user->srec->scan->source.ra_app);break;
      case DecApp   :sscanf(p,"%lf",&user->srec->scan->source.dec_app);break;
      case RaMean   :sscanf(p,"%lf",&user->srec->scan->source.ra_mean);break;
      case DecMean  :sscanf(p,"%lf",&user->srec->scan->source.dec_mean);break;	      case StokesType: sscanf(p,"%hd",&user->stokes_type);break;
      case BurstName:strcpy(user->burst.name,p);break;
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
      case DoFlag   :sscanf(p,"%d",&user->do_flag);break;
      case Thresh   :sscanf(p,"%f",&user->thresh);break;
      case DoBand   :sscanf(p,"%d",&user->do_band);break;
    }
  }

  fclose(fp);
  return 0;
}
/* set up some booking keeping information for each baseline needed at the
   time of copying the data, conversion to FITS etc. Sets up the map between
   the input baseline order and the output random group order.

   jnc mar 2025
   28mar25 fixed bug in loop order

   corrected the offset to reflect that there is no timestamp at the start
   of every record. Instead there is a timestamp at the start of every
   rec_per_slice records. The timestamp is read separately in
   copy_sel_chans() and so does not need to be reflected in the offsets
   here.
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
  // next to each other. 
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
	    visinfo[bs1].flip=0;
	    visinfo[bs1].ant0=ant0;visinfo[bs1].ant1=ant1;
	    visinfo[bs1].band0=b0;visinfo[bs1].band1=b1;
	    visinfo[bs1].drop=0; 
	    visinfo[bs1].off=bs*recl;
	    bs1++;break;
	  }else{
	    if(ant0==a1 && ant1==a0){
	      visinfo[bs1].flip=1;
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
  get the antenna configuration, i.e. names, positions etc.
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
 raw data FITS conversion
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
  corrpar->statime=2.0*corrpar->channels*corrpar->sta/corrpar->clock;
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
  daspar->dut1=37; // seconds as on 16/Mar/25 (WILL BE OVERWRITTEN?)
  
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
  initialize the parameters with sensible values. Almost all of these will be
  overwritten by the actual parameters used at run time. However it is helpful
  to initialize to sensible values since this allows for the code to be
  debugged without having any input files.

  jnc march 2025
*/
int init_user(SvSelectionType *user, char *uparfile, char *antfile){
  VisParType       *vispar=&user->vispar;
  CorrType         *corr=user->corr;
  CorrParType      *corrpar=&corr->corrpar;
  DasParType       *daspar=&corr->daspar;
  RecFileParType   *rfile=&user->recfile;
  ScanRecType      *srec=user->srec;
  ScanInfoType     *scan=srec->scan;
  SourceParType    *source=&scan->source;
  BurstParType     *burst=&user->burst;
  int               i;

  //open the log file
  if((user->lfp=fopen("svfits.log","w"))==NULL)
  { fprintf(stderr,"Cannot open svfits.log\n"); return -1;}

  //set timestamp to the middle of the integration
  user->timestamp_off=0.5*daspar->lta*corrpar->statime;
  user->iatutc = 37 ;  // default
  user->epoch = 2000.0;// default is J2000 coordinates
  user->stokes=2; // RR, LL by default
  /* Set default stokes for RR,LL */
  user->stokes_allow[0] = 0 ;   // RR
  user->stokes_allow[1] = 3 ;   // LL
  user->stokes_allow[2] = user->stokes_allow[3] = -1 ;
  user->stokes_type=0; /* FITS stokes labels are RL */
  user->scans       = 1;
  user->sidebands   = 1;
  user->all_chan    = 0;//ignore burst parms,copy all channels
  user->start_chan  = 0;
  user->chan_inc    = 1;
  user->force_app   = 1; // force coordinates to be treated as apparent
  user->fake_data   = 0; // FOR DEBUGGING!!!! SET TO 0 FOR ACTUAL USE
  user->do_log=20;
  user->do_flag=1;// flag visibilities with outlier amplitudes
  user->update_burst=1;//recompute parms with some as input (see update_burst())
  user->thresh=5.0;// threshold for MAD flagging
  user->do_band=0; // do not do amplitude bandpass
  strcpy(user->fitsfile,"TEST.FITS") ;
  /* default parameters, reset by reading input scanfile*/
  user->hdr->scans=user->scans;
  if((user->lfp=fopen("svfits.log","w"))==NULL)
    {fprintf(stderr,"Unable to open svfits.log\n"); return -1;}
  // initialize the correlator settings
  init_corr(user,antfile); // hardcoded for now
  user->channels=1; //only one output channel by default
  user->antmask=1073741823;//30 antennas (C07 and S05 dropped)
  srec->corr=user->corr;
  srec->scannum=0;
  srec->scan_id=1;
  srec->source_id=1;
  srec->freq_id=1;
  strcpy(scan->proj.code,"TEST"); // replace with GTAC project code
  scan->proj.antmask=daspar->antmask;
  scan->proj.bandmask=daspar->bandmask;
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
  /* default parameteters of the raw data files */
  rfile->nfiles=16;
  strcpy(rfile->path,".");
  for(i=0;i<rfile->nfiles;i++)
    sprintf(rfile->fname[i],"%d.dat",i);// NEED TO UPDATE TO CORRECT CONVENTION
  rfile->rec_per_slice=50; // 50 lta records per slice
  // time duration (sec) of each slice
  rfile->t_slice=rfile->rec_per_slice*daspar->lta*corrpar->statime;
  // time interval (sec) between two successive slices
  rfile->slice_interval=rfile->nfiles*rfile->t_slice;
  rfile->t_start[0]=0;
  for(i=1;i<rfile->nfiles;i++) // start time for data in file
    rfile->t_start[i]=rfile->t_start[i-1]+rfile->t_slice;

  // over ride defaults by parameters given by the user
  if(svuserInp(uparfile,user)) return -1;
  if(user->all_chan) user->channels=corr->daspar.channels;
  //update burst parameters
  if(user->update_burst)update_burst(user);
  // setup the visibility meta data
  if((user->baselines=init_vispar(user))<=0)
    {fprintf(stderr,"No baselines selected!\n"); return -1;}

  return 0;
}
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
int get_rec_num(RecFileParType *rfile, int idx, BurstParType *burst,
	       CorrType *corr, ScanInfoType *scan,int *start_rec,
		int *n_rec,int fake_data){
  double         fb=burst->f/1.0e9;//GHz
  double         burst_width=burst->width;
  double         slice_interval=rfile->slice_interval;
  double         t_slice=rfile->t_slice;
  int            rec_per_slice=rfile->rec_per_slice;
  double         integ=corr->daspar.lta*corr->corrpar.statime;
  double         freq=scan->source.freq[0]/1.0e9; //GHz (first chan)
  double         cw=scan->source.ch_width/1.0e9;//GHz (positive)
  struct timeval tv;
  double         t_burst,t_start,t_rec,fh,fl;
  int            start_slice,rec_num,r;
  FILE          *fp;
  char           fname[LINELEN+PATHLEN+1];
  time_t         unix_time;
  
  if(scan->source.net_sign[0]>0){fh=freq+corr->daspar.channels*cw;fl=freq;}
  else{fh=freq;fl=freq-corr->daspar.channels*cw;}
  
  // time that the trailing end of the burst arrives at the highest freq
  t_burst=burst->t-burst->int_wd/2.0+K0*burst->DM*(1.0/(fh*fh)-1.0/(fb*fb));

  if(!fake_data){
    // get the start time of the file
    strcpy(fname,"");
    if(rfile->path !=NULL){strcpy(fname,rfile->path);strcat(fname,"/");}
    strcat(fname,rfile->fname[idx]);
    fname[strcspn(fname,"\n")]='\0';//remove trailing newline,if any
    if((fp=fopen(fname,"rb"))==NULL)
      {fprintf(stderr,"Unable to open %s\n",fname); return -1;}
    if(fread(&tv,sizeof(struct timeval),1,fp) != 1)
      { fprintf(stderr,"Cannot read the starting time in file %s!\n",
		rfile->fname[idx]);
	return -1;}
    unix_time=tv.tv_sec; // seconds;
    t_start=unix_time_to_ist(unix_time,tv.tv_usec/1.0e6);// IST
    t_start+=idx*rec_per_slice*integ;//timeval is start only for file[0]!
    rfile->t_start[idx]=t_start;
  }else // use the hardcoded start time (debug use only)
    { t_start=rfile->t_start[idx];}
  
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
  *n_rec=r;

  if(!fake_data) fclose(fp);
  


  return 0;
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

   jnc 15/mar/25
*/
int get_nchan(BurstParType *burst,CorrType *corr, ScanInfoType *scan){
  double DM=burst->DM;
  double fb=burst->f/1.0e9; //GHz
  double tb=0; // burst time set to 0, ok since we are just counting channels
  double wd=burst->int_wd;
  double freq=scan->source.freq[0]/1.0e9;//GHz (first chan)
  double cw=scan->source.ch_width/1.0e9; //GHz (positive)
  double integ=corr->daspar.lta*corr->corrpar.statime;
  int    cs,ce,nc,nrec;
  double tr1,tr2,fs,fe,fn,fh;
  
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
  lta record.

  jnc 15/mar/25
*/
int get_chan_num(double trec,BurstParType *burst,CorrType *corr,
		 ScanInfoType *scan,int *cs,int *ce){
  double DM=burst->DM;
  double tb=burst->t;
  double fb=burst->f/1.0e9; //GHz
  double wd=burst->int_wd;
  double freq=scan->source.freq[0]/1.0e9;//GHz (first chan)
  double cw=scan->source.ch_width/1.0e9; //GHz (positive)
  double integ=corr->daspar.lta*corr->corrpar.statime;
  int    c0,c1,c;
  double tr1,tr2,fs,fe,fn,fh,df;
  
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
  // otherwise not reject all channels
  if(c0>0) *cs=c0;else *cs=0;
  if(c1<corr->daspar.channels) *ce=c1;
  else *ce=corr->daspar.channels;
  if(c1<0) // reject all channels
    *cs=*ce=corr->daspar.channels;
  return 0;
}
/*
  compute the uvw coordinates. Modified from gvgetUvw in gvfits

  jnc mar25
*/
void svgetUvw(double tm, double mjd_ref,SourceParType *source, BaseUvwType *uvw)
{ double lst,dec,ha;
  double ch,sh, cd,sd ;
  double bxch,bysh; 
  static double C = 2.99792458e8 ; /* velocity of light in  m/s */
  int k ;

  lst = lmst(mjd_ref + tm/86400) ;
  ha = lst - source->ra_app;
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

  added the option of reading the input file in max_bufsize sized chunks.
  also updated the raw file format, the timestamp is included only once
  at the start of a sice (i.e. at the start of every rec_per_slice records),
  and not for every record.

  jnc apr 2025
*/
int copy_sel_chans(SvSelectionType *user, int idx, char **outbuf,
		   int *restart_rec){
  UvwParType      *uvwpar;
  RecFileParType  *rfile=&user->recfile;
  int              n_slice,rec_per_slice=rfile->rec_per_slice;
  double           slice_interval=rfile->slice_interval;
  CorrType        *corr=user->corr;
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  VisParType      *vispar=&user->vispar;
  VisInfoType     *visinfo=vispar->visinfo;
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  double           epoch,epoch1=user->epoch;
  int              group_size,start_rec,n_rec,n_chan,copied;
  int              rec0,rec1,n_rec1;
  int              r,b,b1,i,date1,c,c0,c1,in_stride;
  double           tm,JD,date2,freq,noise;
  double           integ=corr->daspar.lta*corr->corrpar.statime;
  unsigned long    off,bufsize,max_bufsize,recl,timesize=sizeof(struct timeval);
  BaseUvwType      uvw[MAX_ANTS];
  char            *rbuf,*obuf;
  unsigned short  *in;
  Cmplx3Type      *out;
  struct timeval   tv;
  FILE            *lfp=user->lfp,*fp;
  char             fname[LINELEN+PATHLEN+1];
  time_t           unix_time;
  
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}
  if(idx<0 || idx>rfile->nfiles){
    fprintf(stderr,"Illegal raw file number %d [legal %d - %d]\n",
	    idx,0,rfile->nfiles);
    return -1;
  }
  strcpy(fname,"");
  if(rfile->path !=NULL){strcpy(fname,rfile->path);strcat(fname,"/");}
  strcat(fname,rfile->fname[idx]);
  fname[strcspn(fname,"\n")]='\0';//remove trailing newline,if any
  if((fp=fopen(fname,"rb"))==NULL){
    fprintf(stderr,"Unable to open %s\n",fname); return -1;}

  for (i=0; i<MAX_ANTS; i++)
  { uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }

  // allocate memory for read (timeval comes only once every rec_per_slice
  // records and is read separately
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  if((rbuf=malloc(recl))==NULL){fprintf(stderr,"Malloc Error!\n");return -1;}
  // records containing data - last arg is for fake data
  get_rec_num(rfile,idx,burst,corr,scan,&start_rec,&n_rec,0);
  // max chans containing data - add a safety margin for malloc
  n_chan=(int)rint(1.1*get_nchan(burst,corr,scan)); 
  // each random group contains data for all stokes and one channel
  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes;
  max_bufsize=1024*1024*1024;// copy a maximum of 1 Gigabyte at a time
  bufsize=n_rec*(baselines/stokes)*n_chan*group_size;
  if (bufsize>max_bufsize || restart_rec>0){//process in chunks
    n_rec1=floor((1.0*max_bufsize)/(1.0*bufsize)*n_rec);
    bufsize=n_rec1*(baselines/stokes)*n_chan*group_size;
    rec0=*restart_rec;  rec1=rec0+n_rec1;
    if(rec1>=n_rec){rec1=n_rec;*restart_rec=0;}//last chunk in file
    else{*restart_rec=rec1;}//start of next chunk
  }else {rec0=0;rec1=n_rec;}// process file in one chunk

  obuf=(*outbuf=(char*)malloc(bufsize));//mem returned to calling prog
  if(obuf==NULL) {fprintf(stderr,"Malloc Error\n");return -1;}
  //get the time for the first record of this chunk
  n_slice=floor((start_rec+rec0)/rec_per_slice);
  off=n_slice*sizeof(struct timeval)+n_slice*rec_per_slice*recl;
  if((fseek(fp,off,SEEK_SET)<0)||fread(&tv,timesize,1,fp)!=1){
    fprintf(stderr,"Error reading timestamp from %s\n",rfile->fname[idx]);
    fclose(fp);
    return -1;
  }
  unix_time=tv.tv_sec;
  //middle of record
  tm=unix_time_to_ist(unix_time,tv.tv_usec/1.0e6)+user->timestamp_off;
  tm+=(start_rec+rec0-rec_per_slice*n_slice)*integ;//start time of rec0
  tm+=idx*rec_per_slice*integ;// timeval is start only for file[0]!
  if(user->do_log)
    fprintf(lfp,"File %d  CopyStart=%12.6f Copy Rec %d - %d [Chunk %d - %d]\n",
	    idx,tm,start_rec,start_rec+n_rec,rec0,rec1);

  // convert selected data into random groups
  off=n_slice*timesize+(start_rec+rec0)*recl;//start point of this chunk
  if(fseek(fp,off,SEEK_SET)<0){
    fprintf(user->lfp,"No data selected from File %s\n",rfile->fname[idx]);
    fclose(fp);
    return -1;
  }
  for(copied=0,r=rec0;r<rec1;r++){
    if((start_rec+r)%rec_per_slice ==0){//start of new slice, get timeval
      if(fread(&tv,timesize,1,fp)!=1)
	{fprintf(stderr,"Reached EoF for %s\n",rfile->fname[idx]);return -1;}
      unix_time=tv.tv_sec;
      tm=unix_time_to_ist(unix_time,tv.tv_usec/1.0e6)+user->timestamp_off;
      tm+=idx*rec_per_slice*integ;// timeval is start only for file[0]!
    }else{ tm+=integ;}
    if(fread(rbuf,recl,1,fp)!=1) //get the data
      {fprintf(user->lfp,"Reached EoF %s\n",rfile->fname[idx]); return 0;}
    get_chan_num(tm,burst,corr,scan,&c0,&c1);
    if(user->do_log>10)
      fprintf(lfp,"File %d Rec %d Slice %d Lta %d Time %12.6f Chans %d - %d\n",
	      idx,start_rec+r,(start_rec+r)/50, (start_rec+r)%50, tm,c0,c1);
    if(c0<0 || c0>=channels) continue; // record does not contain burst
    if(c1-c0>n_chan){
      fprintf(user->lfp,"Dropped %d channels in rec %d\n",c1-c0-n_chan,r);
      c1=c0+n_chan;
    }
    svgetUvw(tm,corr->daspar.mjd_ref,source,uvw);
    JD = corr->daspar.mjd_ref + 2400000.5 + tm/86400 ;
    date1 = (int) JD ;
    date2 = JD - date1+user->iatutc/86400 ; //CHECK IATUTC
    for(b=0;b<baselines;b+=vispar->vis_sets){//vis_sets==2
      VisInfoType *vinfo=vispar->visinfo;
      int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
      for(c=c0;c<c1;c++){ // limits for this timestamp
	if(c>=channels)continue;
	freq=freq0+c*df; // frequency of this channel
	off=copied*group_size;
	uvwpar=(UvwParType*)(obuf+off);//uvw for this group
	out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType));//data for this group
	for(b1=b;b1<b+vispar->vis_sets;b1++){
	  //input data for this chan
	  in=(unsigned short*)(rbuf+vinfo[b1].off)+2*c;
	  out->r=half_to_float(in[0]); out->i=half_to_float(in[1]);
	  out->wt=1.0;
	  if(visinfo[b1].flip)out->i=-out->i;
	  if(visinfo[b1].drop)out->wt=-1.0;
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
	epoch=2000.0 + (corr->daspar.mjd_ref - 51544.5)/365.25 ;
	//rotate from date (epoch) coordinates to standard (epoch1==J2000 by
	//default)
	if(epoch1<0.0)epoch1=epoch;
	prenut(uvwpar,corr->daspar.mjd_ref,source->ra_app,source->dec_app,
	       epoch1);
	uvwpar->date1 = date1 ;
	//offset date2 to allow for plotting of the visibility data. The date
	//parameter here is in any case not appropriate for for imaging. Only
	//the u,v,w parameters are correct.
	uvwpar->date2 = date2+(c*integ)/channels ;
	uvwpar->baseline = ant0*256 + ant1 + 257 ;
	uvwpar->su = 1; // only one source in file
	uvwpar->fq = 1 ; // only one frequency id in file
	copied++;
      }
    }
  }

  free(rbuf); // free the read buf; obuf will be freed by calling program
  fclose(fp);
  return copied; // total number of visibility records 
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
  svgetUvw(tm,corr->daspar.mjd_ref,source,uvw); //units metres/C 

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
  like copy_sel_chans() except that no visibility files are read, instead
  some simulated data is generated for debugging.
  
  jnc mar 25
*/
int fake_sel_chans(SvSelectionType *user, int idx, char **outbuf,
		   int *restart_rec){
  UvwParType      *uvwpar;
  RecFileParType  *rfile=&user->recfile;
  CorrType        *corr=user->corr;
  BurstParType    *burst=&user->burst;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  VisParType      *vispar=&user->vispar;
  VisInfoType     *visinfo=vispar->visinfo;
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  double           epoch,epoch1=user->epoch;
  int              recl,group_size,start_rec,n_rec,n_chan,copied;
  int              r,b,b1,i,date1,c,c0,c1,in_stride;
  int              rec0,rec1,n_rec1;
  double           tm,tm_m,JD,date2,freq,noise;
  double           integ=corr->daspar.lta*corr->corrpar.statime;
  unsigned long    off,bufsize,max_bufsize;
  BaseUvwType      uvw[MAX_ANTS];
  char            *rbuf,*obuf;
  unsigned short        *in;
  Cmplx3Type      *out;
  FILE            *lfp=user->lfp;
  
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}

  for (i=0; i<MAX_ANTS; i++)
  { uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }

  if(idx<0 || idx>rfile->nfiles){
    fprintf(stderr,"Illegal raw file number %d [legal %d - %d]\n",
	    idx,0,rfile->nfiles);
    return -1;
  }
  // allocate memory for read
  recl=corr->daspar.baselines*corr->daspar.channels*sizeof(float);
  if((rbuf=malloc(recl))==NULL){fprintf(stderr,"Malloc Error!\n");return -1;}
  
  // get recnums without actuall reading the data (last arg is 1 for fakedata)
  get_rec_num(rfile,idx,burst,corr,scan,&start_rec,&n_rec,1);
  // add a safety margin for malloc
  n_chan=(int)rint(1.1*get_nchan(burst,corr,scan)); 
  // each random group contains data for all stokes and one channel
  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes;
  max_bufsize=1024*1024*1024;// copy a maximum of a Gigabyte at a time
  bufsize=n_rec*(baselines/stokes)*n_chan*group_size;
  if (bufsize>max_bufsize || restart_rec>0){//process in chunks
    n_rec1=floor((1.0*max_bufsize)/(1.0*bufsize)*n_rec);
    bufsize=n_rec*(baselines/stokes)*n_chan*group_size;
    rec0=*restart_rec;
    rec1=rec0+n_rec1;
    if(rec1>=n_rec){
      rec1=n_rec;// last chunk of file
      *restart_rec=0;
    }else{
      *restart_rec=rec1;//start point for next chunk to copy
    }
  }else{
    rec0=0;
    rec1=n_rec;
  }
  *outbuf=(char*)malloc(n_rec*(baselines/stokes)*n_chan*group_size);
  obuf=*outbuf;
  if(obuf==NULL) {fprintf(stderr,"Malloc Error\n");return -1;}
  if(user->do_log)
    fprintf(lfp,"File %d  Copy Records %d - %d [chunk: %d-%d]\n",
	    idx,start_rec,start_rec+n_rec,start_rec+rec0,start_rec+rec1);
  tm=rfile->b_start[idx];
  for(copied=0,r=rec0;r<rec1;r++){
    get_chan_num(tm,burst,corr,scan,&c0,&c1);
    if(user->do_log>10)
      fprintf(lfp,"File %d Rec %d Slice %d Lta %d Time %12.6f Copy Chans %d - %d\n",idx,start_rec+r,(start_rec+r)/50, (start_rec+r)%50, tm,c0,c1);
    if(c0<0 || c0>=channels) continue; // record does not contain burst
    if(c1-c0>n_chan){
      fprintf(user->lfp,"Dropped %d channels in rec %d\n",c1-c0-n_chan,r);
      c1=c0+n_chan;
    }
    tm_m=tm+user->timestamp_off;// middle of rec
    if(simulate_visibilities(user, tm_m,rbuf,c0,c1))return -1;
    svgetUvw(tm_m,corr->daspar.mjd_ref,source,uvw);
    JD = corr->daspar.mjd_ref + 2400000.5 + tm_m/86400 ;
    date1 = (int) JD ;
    date2 = JD - date1+user->iatutc/86400 ; //CHECK IATUTC
    for(b=0;b<baselines;b+=vispar->vis_sets){//vis_sets==2
      VisInfoType *vinfo=vispar->visinfo;
      int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
      for(c=c0;c<c1;c++){ // limits for this timestamp
	if(c>=channels)continue;
	freq=freq0+c*df; // frequency of this channel
	off=copied*group_size;
	uvwpar=(UvwParType*)(obuf+off);//uvw for this group
	out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType));//data for this group
	for(b1=b;b1<b+vispar->vis_sets;b1++){
	  //input data for this chan
	  in=(unsigned short*)(rbuf+vinfo[b1].off)+2*c;
	  out->r=half_to_float(in[0]);
	  out->i=half_to_float(in[1]);
	  noise=drand48()-0.5;out->r +=0.1*noise;
	  noise=drand48()-0.5;out->i +=0.1*noise;
	  if(visinfo[b1].flip)out->i=-out->i;
	  if(visinfo[b1].drop)out->wt=-1.0;else out->wt=1.0;
	  out++;//now copy same chan for next stokes
	 }
	//NEED TO compute uvw and ROTATE TO J2000!!!!
	// scale uv co-ordinates so that the values come out right inside of
	// AIPS/CASA. In the imaging programs the number here is multiplied
	// by the frequency of the channel inorder to get the uv coordinate
	// for that channel). In the FITS file, all of the data is put in
	// the first channel, i.e. at frequency freq0. The scaling here
	// corrects for this.
	uvwpar->u = (uvw[ant1].u-uvw[ant0].u)*(freq/freq0);
	uvwpar->v = (uvw[ant1].v-uvw[ant0].v)*(freq/freq0);
	uvwpar->w = (uvw[ant1].w-uvw[ant0].w)*(freq/freq0);
	epoch=2000.0 + (corr->daspar.mjd_ref - 51544.5)/365.25 ;
	//rotate from date (epoch) coordinates to standard (epoch1==J2000 by
	//default)
	if(epoch1<0.0)epoch1=epoch;
	prenut(uvwpar,corr->daspar.mjd_ref,source->ra_app,source->dec_app,
	       epoch1);
	uvwpar->date1 = date1 ;
	uvwpar->date2 = date2 ;
	uvwpar->baseline = ant0*256 + ant1 + 257 ;
	uvwpar->su = 1; // only one source in file
	uvwpar->fq = 1 ; // only one frequency id in file
	copied++;
      }
    }
    if((start_rec+r)%rfile->rec_per_slice==rfile->rec_per_slice-1)
      tm+=rfile->slice_interval-rfile->t_slice;
    else tm+=integ;
  }

  free(rbuf); // free the read buf; obuf will be freed by calling program
  return copied; // total number of visibility records 
}
int clip(char *visbuf, SvSelectionType *user, int groups){
  double         thresh=user->thresh;
  int            stokes=user->stokes;
  Cmplx3Type    *vis;
  long           off=0;
  float         *amp,*med_a,*mad_a,med,mad,a,dum;
  int            sets,set_size,sel_sets,group_size;
  int            i,j,k,flagged=0;

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

  //flag the data
  for(i=0;i<groups;i++){
    off=i*group_size+sizeof(UvwParType);
    vis=(Cmplx3Type*)(visbuf+off);
    a=sqrt(vis[0].r*vis[0].r+vis[0].i*vis[0].i);
    if(fabs(a-med)>thresh*mad) {vis[0].wt=-1;flagged++;}
    a=sqrt(vis[1].r*vis[1].r+vis[1].i*vis[1].i);
    if(fabs(a-med)>thresh*mad){vis[1].wt=-1;flagged++;}
  }

  free(amp);free(med_a);free(mad_a);

  return flagged;
}
/*
  function to compute the bandpass. In detail, it returns in the user->bpass
  structure:(1) the average (over all records in the slice - rbuf is expected to
  point to a buffer holding all the records in a given slice) of the real
  and imaginary part separately for each channel of each baselines selected
  for output. This average (arithmentic mean) can be subtracted from each
  visibility - the idea being that this would get rid of some RFI at least.
  (Also gets rid of some signal - but the intention is that when do something
  similar for the burst (in another as yet [13/apr/25] to be written function)
  we will exclude visibilities with signal from the average. (2) the average
  bandpass amplitude. It is ensured that the average amplitude bandpass has
  some sensible value in it, even if some channels have no data. The amplitude
  bandpass is normalized using the median value across all channels.

  This particular function is mainly meant for debugging use, for e.g. to check
  that the calibration is working properly, etc.

  jnc apr 2025
*/
  
int make_bpass(SvSelectionType *user, char *rbuf){
  int              rec_per_slice=user->recfile.rec_per_slice;
  CorrType        *corr=user->corr;
  VisParType      *vispar=&user->vispar;
  VisInfoType     *visinfo=vispar->visinfo;
  BpassType       *bpass=&user->bpass;
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  unsigned long    off,recl;
  int              r,n,b,c;
  unsigned short  *in;
  Cmplx3Type      *mean;
  float            re,im,median,mad,*amp;

  // allocate memory for the mean over records as well as the normalized
  // amplitude. Calling program should free memory when done
  for(b=0;b<baselines;b++){
    if((bpass->mean[b]=(Cmplx3Type*)malloc(channels*sizeof(Cmplx3Type)))==NULL)
      {fprintf(stderr,"Malloc error in make_bpass\n"); return -1;}
    if((bpass->abp[b]=(float*)malloc(channels*sizeof(float)))==NULL)
      {fprintf(stderr,"Malloc error in make_bpass\n"); return -1;}
  }

  // compute the mean, normalized amplitude bandpass for selected channels
  for(b=0;b<baselines;b++){
    mean=bpass->mean[b];
    for(c=0;c<channels;c++){
      mean[c].r=mean[c].i=0.0;
      for(n=0,r=0;r<rec_per_slice;r++){
	off=r*recl+visinfo[b].off;
	in=(unsigned short*)(rbuf+off)+2*c;
	re=half_to_float(in[0]);
	im=half_to_float(in[1]);
	if(isfinite(re) && isfinite(im))
	  {mean[c].r+=re; mean[c].i+=im;n++;}
      }
      if(n>0)
	{mean[c].r=mean[c].r/n;mean[c].i=mean[c].i/n;mean[c].wt=1.0;}
      else mean[c].wt=-1;
    }
    amp=bpass->abp[b];
    for(c=0;c<channels;c++)
      if(mean[c].wt>0)
	amp[c]=sqrt(mean[c].r*mean[c].r+mean[c].i*mean[c].i);
      else
	amp[c]=-1.0;
    //interpolate over flagged channels
    if(amp[0]<0.0){
      for(c=1;c<channels;c++)
	if(amp[c]>0.0)
	  {amp[0]=amp[c];break;}
    }
    for(c=1;c<channels;c++)
      if(amp[c]<0.0)amp[c]=amp[c-1];
    if(amp[0]<0.0)
      {fprintf(stderr,"All channels flagged, cannot normalize!\n"); return -1;}
    robust_stats(channels,amp,&median,&mad);
    for(c=0;c<channels;c++)
      amp[c]=amp[c]/median;
  }
  return 0;
}
/*
  function to copy all channels in a given slice. This then creates a "normal"
  UV file, i.e. there is no scaling of UV coordinates etc. Mainly for debugging
  use, i.e. for e.g. for checking that the bandpass calibration is working,
  etc.

  jnc apr 2025
*/
int copy_all_chans(SvSelectionType *user, int idx, char **outbuf,
		   int restart_slice){
  UvwParType      *uvwpar;
  RecFileParType  *rfile=&user->recfile;
  int              n_slice,rec_per_slice=rfile->rec_per_slice;
  double           slice_interval=rfile->slice_interval;
  CorrType        *corr=user->corr;
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  double           freq0=source->freq[0];
  double           df=source->ch_width*source->net_sign[0];//signed
  VisParType      *vispar=&user->vispar;
  VisInfoType     *visinfo=vispar->visinfo;
  int              channels=corr->daspar.channels;//input channels
  int              baselines=user->baselines;//output baselines
  int              stokes=user->stokes;
  double           epoch,epoch1=user->epoch;
  int              group_size,start_rec,n_rec,copied;
  int              rec0,rec1,n_rec1;
  int              r,b,b1,i,date1,c;
  double           tm,JD,date2,freq;
  double           integ=corr->daspar.lta*corr->corrpar.statime;
  unsigned long    off,bufsize,max_bufsize,recl,timesize=sizeof(struct timeval);
  BaseUvwType      uvw[MAX_ANTS];
  char            *rbuf,*obuf;
  unsigned short  *in;
  Cmplx3Type      *out,*mean,*mean0;
  float           *abp,*abp0;
  struct timeval   tv;
  FILE            *lfp=user->lfp,*fp;
  char             fname[LINELEN+PATHLEN+1];
  time_t           unix_time;
  
  if (baselines < 1)
    {fprintf(stderr,"Illegal Nbaselines %d\n",baselines); return -1 ;}
  if ((stokes != 1) && (stokes != 2) &&(stokes !=4))
    {fprintf(stderr,"Illegal Nstokes %d (Legal 1/2/4)\n",stokes);return -1;}
  if(idx<0 || idx>rfile->nfiles){
    fprintf(stderr,"Illegal raw file number %d [legal %d - %d]\n",
	    idx,0,rfile->nfiles);
    return -1;
  }
  strcpy(fname,"");
  if(rfile->path !=NULL){strcpy(fname,rfile->path);strcat(fname,"/");}
  strcat(fname,rfile->fname[idx]);
  fname[strcspn(fname,"\n")]='\0';//remove trailing newline,if any
  if((fp=fopen(fname,"rb"))==NULL){
    fprintf(stderr,"Unable to open %s\n",fname); return -1;}

  for (i=0; i<MAX_ANTS; i++)
  { uvw[i].bx = user->corr->antenna[i].bx ;
    uvw[i].by = user->corr->antenna[i].by ;
    uvw[i].bz = user->corr->antenna[i].bz ;
  }

  // allocate memory for read (timeval comes only once every rec_per_slice
  // records and is read separately
  recl=corr->daspar.baselines*channels*sizeof(float);
  if((rbuf=malloc(rec_per_slice*recl))==NULL){
    fprintf(stderr,"Malloc Error allocated %ld bytes!\n",rec_per_slice*recl);
      return -1;
  }
  group_size=sizeof(UvwParType)+sizeof(Cmplx3Type)*stokes*channels;
  bufsize=rec_per_slice*(baselines/stokes)*group_size;
  obuf=(*outbuf=(char*)malloc(bufsize));//mem returned to calling prog
  if(obuf==NULL) 
    {fprintf(stderr,"Malloc Error allocating %ld bytes\n",bufsize);return -1;}
  //get the time for the first record of this chunk
  off=restart_slice*sizeof(struct timeval)+restart_slice*rec_per_slice*recl;
  if((fseek(fp,off,SEEK_SET)<0)||fread(&tv,timesize,1,fp)!=1){
    fprintf(stderr,"Error reading timestamp from %s\n",rfile->fname[idx]);
    fclose(fp);
    return -1;
  }
  unix_time=tv.tv_sec;
  //middle of record
  tm=unix_time_to_ist(unix_time,tv.tv_usec/1.0e6)+user->timestamp_off;
  tm+=idx*rec_per_slice*integ;// timeval is start only for file[0]!
  if(user->do_log)
    fprintf(lfp,"File %d  CopyStart=%12.6f Copy Rec %d - %d [Chunk %d - %d]\n",
	    idx,tm,start_rec,start_rec+n_rec,rec0,rec1);
  if(fread(rbuf,recl,rec_per_slice,fp)!=rec_per_slice) //get the data
    {fprintf(user->lfp,"Reached EoF %s\n",rfile->fname[idx]); return 0;}
  if(user->do_band)    //compute the bandpass
    {if(make_bpass(user,rbuf)) return -1;}
  else{// calibration arrays that won't affect data values
    if((mean0=(Cmplx3Type*)malloc(channels*sizeof(Cmplx3Type)))==NULL)
      {fprintf(stderr,"Malloc error\n");return -1;}
    if((abp0=(float*)malloc(channels*sizeof(float)))==NULL)
      {fprintf(stderr,"Malloc error\n");return -1;}
    for(c=0;c<channels;c++)
      {mean0[c].r=mean0[c].i=0.0;abp0[c]=mean0[c].wt=1.0;}
  }
	    
  // convert selected data into random groups
  for(copied=0,r=0;r<rec_per_slice;r++){
    svgetUvw(tm,corr->daspar.mjd_ref,source,uvw);
    JD = corr->daspar.mjd_ref + 2400000.5 + tm/86400 ;
    date1 = (int) JD ;
    date2 = JD - date1+user->iatutc/86400 ; //CHECK IATUTC
    for(b=0;b<baselines;b+=vispar->vis_sets){//vis_sets==2
      VisInfoType *vinfo=vispar->visinfo;
      int ant0=vinfo[b].ant0,ant1=vinfo[b].ant1; //same for all base in set
      off=copied*group_size;
      uvwpar=(UvwParType*)(obuf+off);//uvw for this group
      out=(Cmplx3Type*)(obuf+off+sizeof(UvwParType));//data for this group
      for(b1=b;b1<b+vispar->vis_sets;b1++){
	off=r*recl+vinfo[b1].off;
	if(user->do_band){//actual calibration
	  mean=user->bpass.mean[b1]; abp=user->bpass.abp[b1];
	}else{ mean=mean0; abp=abp0;}//dummy values
	in=(unsigned short*)(rbuf+off);
	for(c=0;c<channels;c++){ 
	  //input data for this chan
	  out->r=half_to_float(in[2*c]);
	  out->i=half_to_float(in[2*c+1]);
	  out->r=(out->r-mean[c].r)/abp[c];
	  out->i=(out->i-mean[c].i)/abp[c];
	  out->wt=1.0;
	  if(visinfo[b1].flip)out->i=-out->i;
	  if(visinfo[b1].drop)out->wt=-1.0;
	  if(!(isfinite(out->r) &&isfinite(out->i)))
	    out->wt=-1.0;
	  out++;//now copy same chan for next stokes
	}
      }
      uvwpar->u = (uvw[ant1].u-uvw[ant0].u);
      uvwpar->v = (uvw[ant1].v-uvw[ant0].v);
      uvwpar->w = (uvw[ant1].w-uvw[ant0].w);
      epoch=2000.0 + (corr->daspar.mjd_ref - 51544.5)/365.25 ;
      //rotate from date (epoch) coordinates to standard (epoch1==J2000 by
      //default)
      if(epoch1<0.0)epoch1=epoch;
      prenut(uvwpar,corr->daspar.mjd_ref,source->ra_app,source->dec_app,
	     epoch1);
      uvwpar->date1 = date1 ;
      uvwpar->date2 = date2 ;
      uvwpar->baseline = ant0*256 + ant1 + 257 ;
      uvwpar->su = 1; // only one source in file
      uvwpar->fq = 1 ; // only one frequency id in file
      copied++;
      tm=tm+integ;
    }
  }

  free(rbuf); // free the read buf; obuf will be freed by calling program
  fclose(fp);
  if(user->do_band){
    for(b=0;b<baselines;b++)
      {free(user->bpass.mean[b]); free(user->bpass.abp[b]);}
  }else{
    free(mean0); free(abp0);
  }
  return copied; // total number of visibility records 
}
