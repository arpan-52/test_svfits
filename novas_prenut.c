#include"novas.h"
#include"novascon.h"
#include"solarsystem.h"
#include"nutation.h"

#include "svio.h"

/*
  function to compute the J2000 from the mean coordinates using the novas
   routines
*/
int novas_mean2j2000(double ram,double decm,double mjd,double iatutc, 
	       double *ra,double *dec){
  double jd_tdb1,jd_tdb2;
  double pos1[3],pos2[3],dist;

  dist=1.0e7; // some large distance (units AU)
  // computing TT not TDB (more than good enough)
  jd_tdb1=2400000.5+mjd+(iatutc+32.184)/86400.0; 
  radec2vector(ram*(12.0/M_PI),decm*(180.0/M_PI),dist,pos1);
  jd_tdb2=2451545.0; // J2000
  if(precession(jd_tdb1,pos1,jd_tdb2,pos2)) return -1;
  vector2radec(pos2,ra,dec);
  *ra=(*ra)*(M_PI/12.0); // radians
  *dec=(*dec)*(M_PI/180.0); // radians
  return 0;
}
/*
  Function to compute the J2000 coordinates using the novas routines from
  the apparent coordinates.
  jnc may 2025
*/
void novas_app2j2000(double rap, double decp, double mjd, double iatutc,
	       double *ra,double *dec)
{
  double    jd_tt;
  short int accuracy;

  accuracy=1; // less accurate, fine for GMRT
  //TT = UTC+IATUTC+32.184s
  jd_tt=2400000.5+mjd+(iatutc+32.184)/86400;// now in TT
  rap=rap*(180.0/M_PI)/15.0;//rad->hours
  decp=decp*(180.0/M_PI); // rad->degrees
  mean_star(jd_tt,rap,decp,accuracy,ra,dec); //for app->icrs (j2000)
  *ra=(*ra)*15.0*(M_PI/180.0); // hours->radians
  *dec=(*dec)*(M_PI/180.0); // degrees -> radians

  return;
}  
/*
  Function to rotate the antenna based (uvw) coordinates. Specifically the
  coordinates are rotated twice, the first to undo the nutation correction
  and then to undo the precession. Uses the novasc library calls.

  novasc routines work in the equitorial frame, so the vectors need to be
  rotated to the equitorial frame, corrected for precession and nutation and
  then rotated back to the interferometric (i.e. uvw) frame. See init_mat()
  for the matrix initializations for this. It is assumed that this
  initialization has been done before this routine is called. Initialization
  assumes that there is no change of the phase centre during
  the observation.

  To get to J2000 co-ordinates one also has to correct for the aberration,
  this is treated as a shift and not a rotation, and is corrected by
  puting the correct J2000 coordinates for the reference direction.


  jnc may 2025
*/
#define MAT_VEC(A,v,v1){for(i1=0;i1<3;i1++){for(v1[i1]=0,j1=0;j1<3;j1++){v1[i1]+=A[i1][j1]*v[j1];}}}
int novas_prenut_uvw(SvSelectionType *user,BaseUvwType *uvw, double mjd){
  double jd0,jd1,p0[3],p1[3];
  int    i1,j1,k;// i1,j1 used in MAT_VEC, do *not* reuse.
  short  int direction=1;// undo nutation to give mean coordinates
  short  int accuracy=1; //reduced accuracy (good enough for GMRT)

  jd0=2451545.0; //J2000
  jd1=2400000.5+mjd; 
  for(k=0;k<MAX_ANTS;k++){
    p1[0]=uvw[k].u; p1[1]=uvw[k].v; p1[2]=uvw[k].w;
    MAT_VEC(user->i2eapp,p1,p0);// rotate from UVW to equitorial frame
    nutation(jd1,direction,accuracy,p0,p1);//undo nutation
    precession(jd1,p1,jd0,p0);// undo precession (i.e. now J2000)
    MAT_VEC(user->emean2i,p0,p1);// rotate from equitorial to UVW frame
    uvw[k].u=p1[0]; uvw[k].v=p1[1]; uvw[k].w=p1[2];
  }

  return 0;
}
/*
  rotates the visibility uvw coordinates, unlike prenut_uvw which
  rotates the per antenna uvw coordinates. wrote this as part of
  debugging to see if it makes a difference

  jnc 11/sep/25
*/
int novas_prenut_vis(SvSelectionType *user,UvwParType *uvw, double mjd){
  double jd0,jd1,p0[3],p1[3];
  int    i1,j1,k;// i1,j1 used in MAT_VEC, do *not* reuse.
  short  int direction=1;// undo nutation to give mean coordinates
  short  int accuracy=1; //reduced accuracy (good enough for GMRT)

  jd0=2451545.0; //J2000
  jd1=2400000.5+mjd; 
  p1[0]=uvw->u; p1[1]=uvw->v; p1[2]=uvw->w;
  MAT_VEC(user->i2eapp,p1,p0);// rotate from UVW to equitorial frame
  nutation(jd1,direction,accuracy,p0,p1);//undo nutation
  precession(jd1,p1,jd0,p0);// undo precession (i.e. now J2000)
  MAT_VEC(user->emean2i,p0,p1);// rotate from equitorial to UVW frame
  uvw->u=p1[0]; uvw->v=p1[1]; uvw->w=p1[2];

  return 0;
}
/*
  function to initalize the rotation matrix for rotating the visibilities
  to a new phase centre, as well as the lmn coordinates of the new phase
  centre in the old coordinate system. Used in vis_recentre().

  All uvw calculations are done in J2000, since it is assumed that the input
  uvw coordinates to vis_recentre are in J2000. The actual rotation operation
  is done on the visibilities after conversion to random group format. See
  vis_recentre() for more details.

  The lmn components of the vector connecting the new phase centre with the
  old one are computed in both the J2000 frame, as well as the original
  apparent coordinates frame. The former is used in vis_recentre(), while the
  latter is meant for use in make_postcorr_beam(), where no rotation to J2000
  is done.

  jnc may 2025

  fixed a inversion in the matrix multiplication order

  jnc 14june 2025
*/
#define MAT_MAT(A,B,C){for(i1=0;i1<3;i1++){for(j1=0;j1<3;j1++){for(C[i1][j1]=0,k1=0;k1<3;k1++){C[i1][j1]+=A[i1][k1]*B[k1][j1];}}}}
void init_mat(SvSelectionType *user){
  ScanInfoType    *scan=user->srec->scan;
  SourceParType   *source=&scan->source;
  int              stokes=user->stokes;
  BurstParType    *burst=&user->burst;
  double           r0,d0,r1,d1,t;
  double           m0[3][3],m1[3][3];
  int              i1,j1,k1;
  
  // phase and burst beam centres in J2000 coordinates
  novas_app2j2000(source->ra_app,source->dec_app,user->recfile.mjd_ref,user->iatutc,
	    &r0,&d0); 
  novas_app2j2000(burst->ra_app,burst->dec_app,user->recfile.mjd_ref,user->iatutc,
	    &r1,&d1); 

  //matrix to rotate from interferometric (uvw) axis to date equitorial axis
  /* rotate from source(dec) to equitorial plane*/
  t=(M_PI/2-source->dec_app);
  m0[0][0]=1.0; m0[0][1]=0.0;     m0[0][2]= 0.0;
  m0[1][0]=0.0; m0[1][1]=cos(t);  m0[1][2]=-sin(t);
  m0[2][0]=0.0; m0[2][1]=sin(t);  m0[2][2]= cos(t);
  /* rotate from source(ra) to ra=0.0 */
  t=(M_PI/2.0+source->ra_app);
  m1[0][0]=cos(t); m1[0][1]=-sin(t); m1[0][2]=0.0;
  m1[1][0]=sin(t); m1[1][1]= cos(t); m1[1][2]=0.0;
  m1[2][0]=0.0;    m1[2][1]= 0.0;    m1[2][2]=1.0;
  MAT_MAT(m1,m0,user->i2eapp); // net rotation

  //matrix to rotate from J2000 equitorial axis to interferometric axis
  /* rotate from equitorial plane to source(dec) (NB indices transposed)*/
  t=(M_PI/2-d0);
  m0[0][0]=1.0; m0[1][0]=0.0;     m0[2][0]= 0.0;
  m0[0][1]=0.0; m0[1][1]=cos(t);  m0[2][1]=-sin(t);
  m0[0][2]=0.0; m0[1][2]=sin(t);  m0[2][2]= cos(t);
  /* rotate from ra=0.0 to source(ra) (NB indices transposed)*/
  t=(M_PI/2.0+r0);
  m1[0][0]=cos(t); m1[1][0]=-sin(t); m1[2][0]=0.0;
  m1[0][1]=sin(t); m1[1][1]= cos(t); m1[2][1]=0.0;
  m1[0][2]=0.0;    m1[1][2]= 0.0;    m1[2][2]=1.0;
  MAT_MAT(m0,m1,user->emean2i); //net rotation

  if(!user->recentre) return;
  
  //set up the matrix to rotate visibility to new centre
  user->rmat[0][0]=cos(r1-r0);
  user->rmat[0][1]=sin(d0)*sin(r1-r0);
  user->rmat[0][1]=-cos(d0)*sin(r1-r0);
    
  user->rmat[1][0]=sin(d0)*sin(r1-r0);
  user->rmat[1][1]=sin(d0)*sin(d1)*cos(r1-r0)+cos(d0)*cos(d1);
  user->rmat[1][2]=sin(d0-d1);
    
  user->rmat[2][0]=cos(d1)*sin(r1-r0);
  user->rmat[2][1]=sin(d1)*cos(d0)-cos(d1)*sin(d0)*cos(r1-r0);
  user->rmat[2][2]=cos(d0)*cos(d1)*cos(r1-r0)+sin(d0)*sin(d1);

  //set up the source direction cosines in the original uvw basis at J2000.
  // Used to correct the visibility phase in vis_recentre()
  user->lmn[0]=cos(d1)*sin(r1-r0);
  user->lmn[1]=sin(d1)*cos(d0)-cos(d1)*sin(d0)*cos(r1-r0);
  user->lmn[2]=sqrt(1.0-user->lmn[0]*user->lmn[0]-user->lmn[1]*user->lmn[1]);

  //set up the source direction cosines in the original uvw basis at DATE
  // Used to correct the visibility phase in make_postcorr_beam()
  r0=source->ra_app;d0=source->dec_app;
  r1=burst->ra_app; d1=burst->dec_app;
  user->lmn_a[0]=cos(d1)*sin(r1-r0);
  user->lmn_a[1]=sin(d1)*cos(d0)-cos(d1)*sin(d0)*cos(r1-r0);
  user->lmn_a[2]=sqrt(1.0-user->lmn_a[0]*user->lmn_a[0]-
		      user->lmn_a[1]*user->lmn_a[1]);
  

  fprintf(user->lfp,"Rotating the visibilities to the Burst Beam Centre\n");

  return;
}
