#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "svio.h"

extern void sla_nut_(double *mjd, double *a); 
extern void sla_prec_(double *epoch0,double *epoch1, double *a); 
extern double sla_epj_(double *mjd);
extern void   sla_amp_(double *ra_app, double *dec_app, double *mjd,
		       double *epoch1,double *ra_mean, double *dec_mean);

#define TINY 1.0e-9
/*
  This function rotates the uvw coordinates of the individual antennas, and
  was testing only. Not used.
  JNC 4/Oct/25
*/
void sla_prenut(BaseUvwType *uvw,double mjd,double ra_app,double dec_app,
	       double epoch1)
{ double  epoch,v[3],v1[3],a[9],nm[3][3],pm[3][3],t,ra_mean,dec_mean;
  static double  p[3][3],p1[3][3],p2[3][3],p3[3][3],rm[3][3];
  int    i,j,k,kk;
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
    fprintf(stderr,"EPOCH=%f EPOCH1=%f\n",epoch,epoch1);
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


  for(kk=0;kk<MAX_ANTS;kk++){
    v[0]=uvw[kk].u;v[1]=uvw[kk].v;v[2]=uvw[kk].w;
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
    uvw[kk].u=v1[0];uvw[kk].v=v1[1];uvw[kk].w=v1[2]; 
  }

  return;
}
/*
  This function is identical to the prenut() function in gvfits.

  Not thread safe since it uses static variables, but in the svfits application
  the source and epoch never changes, so the static variable logic is not
  critical.

    JNC 4/Oct/25
*/
void sla_prenut_vis(UvwParType *uvw,double mjd,double ra_app,double dec_app,
		   double epoch1)
{ double  epoch,v[3],v1[3],a[9],nm[3][3],pm[3][3],t,ra_mean,dec_mean;
  static double  p[3][3],p1[3][3],p2[3][3],p3[3][3],rm[3][3];
  int    i,j,k,kk;
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
      fprintf(stderr,"MJD=%f EPOCH=%f EPOCH1=%f\n",mjd,epoch,epoch1);
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


  v[0]=uvw->u;v[1]=uvw->v;v[2]=uvw->w;
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
  uvw->u=v1[0];uvw->v=v1[1];uvw->w=v1[2]; 
  return;
}

#undef TINY

