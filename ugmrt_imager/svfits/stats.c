/*
    Copyright 2013,2014 Jayaram N. Chengalur

    This file is part of flagcal. flagcal is free software: you can 
    redistribute it and/or modify it under the terms of the GNU General 
    Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version.

    flagcal is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
    for more details.

    You should have received a copy of the GNU General Public License
    along with flagcal.  If not, see <http://www.gnu.org/licenses/>.

    Correspondence regarding flagcal should be addressed to
    Jayaram N Chengalur (chengalur@ncra.tifr.res.in)
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void straight_sort(float *x, int n){
  /* traditional sorting routine. This is an n^2 algorithm, but
     is generally faster than quicksort for short arrays
     jnc 06/sep/12
  */
  int   i,j;
  float tmp;

  for(i=1;i<n;i++){
    tmp = x[i];
    for(j=i-1;j>-1 && x[j]>tmp;j--)
      x[j+1] = x[j];
    x[j+1]=tmp;
  }

}
#define QS_SWAP(u,v) {QS_tmp=u;u=v;v=QS_tmp;}
float quickselect(float *x, int n, int k, int *status){
/* C. A. R Hoare's partition based selection algorithm 
   jnc 05/sep/12
*/
     
  int    l,r,m,p,s,i,neq;
  float  val,QS_tmp;
  int    ir=11;   /* seed for random number generator */
  int    niter=0;
  double maxiter;

  *status=0;
  maxiter=pow(n,3);
  
  if(n <= 0){
    *status=-1;
    fprintf(stderr,"quickselect called with n=%d\n!",n);
    return 0.0;
  }
  if(n==1) return x[0];
  if(n==2){
    if(k==1) 
      return x[0]>x[1]? x[1]:x[0];
    else
      return x[0]>x[1]? x[0]:x[1];
  }

  l=0;     /* start position for search */
  r=n-1;   /* end position for search   */
  m=n;     /* length of array to search */
  while(1){
    ir = (ir*211+1663)%7875; /* quick and dirty random num, see NR */
    p = l +m*(ir/7875.0); /* pivot position */
    val = x[p];
    s   = l;
    neq = 0;
    QS_SWAP(x[p],x[r]);
    for(i=l;i<r;i++){
      if(x[i]<=val){
	if(x[i]==val)neq++;
	if(i!=s) QS_SWAP(x[s],x[i]);
	s++;
      }
    }
    if(s==k)return x[r];
    if(s>k){
      if(s-neq <= k) return val;
      r=s;
    }
    if(s<k) l=s;
    m=r-l+1;
    QS_SWAP(x[r],x[s]);
    niter++;
    if(niter > maxiter){
      fprintf(stderr,"Used more than n^3 (n=%d)iterations - quitting now\n",n);
      *status=-1;
      return 0.0;
    }
  }

  /* should never get here! */
  fprintf(stderr,"Reached unexpected branch of code\n");
  *status=-1;

  return 0.0;
}
#undef QS_SWAP
float quick_median(float *x, int n, int *status){
  /* function to return the median. Decides which algorithm
     (straight_sort or quickselect) to use, depending on the
     array length. straight_sort is fastest for short arrays
     while quickselect is better for long arrays. For even
     length arrays one would have to call quickselect twice
     to determine the median, so in this case, the straight_sort 
     is faster for longer arrays than in the case
     of odd length arrays. These two algorithms between them
     seem to be sufficient, and there does not seem to be
     any great advantage in trying also quicksort in intermediate
     length cases.
    
     The exact length at which the switchover from straight
     sort to quickselect needs to be done is probably
     machine dependent. The current values have been
     set on orca.

     jnc 08/sep/12

  */
  float med;

  if(n%2){ /* odd number of elements */
    if(n<36){
      straight_sort(x,n);
      med = x[n/2];
      *status=0;
    }else{
      med=quickselect(x,n,n/2,status);
      if(*status) return 0.0;
    }
  }else{ /* even number of elements */
    if(n<73){
      straight_sort(x,n);
      med = 0.5*(x[n/2]+x[n/2-1]);
      *status=0;
    }else{
      med  = quickselect(x,n,n/2,status);
      if(*status) return 0.0;
      med +=quickselect(x,n,n/2-1,status);
      if(*status) return 0.0;
      med *=0.5;
    }
  }

  return med;
}

int robust_stats(int n, float *x, float *med, float *mad){
  float *x1;
  int    i;
  int    status;


  if(n<0){
    fprintf(stderr,"Negative array length %d illegal\n",n);
    return 1;
  }

  if(n==0){*med=*mad=0.0; return 0;}
  if(n==1){*med=x[0];*mad=0; return 0;}

  if((x1=(float*)malloc(n*sizeof(float)))==NULL)
    {fprintf(stderr,"malloc error\n"); return -1;}

  for(i=0;i<n;i++)x1[i]=x[i];
  *med = quick_median(x1, n, &status);
  if(status) return status;
  for(i=0;i<n;i++)x1[i]=fabs(x[i]-(*med));
  *mad = quick_median(x1, n, &status);
  if(status)return status;

  free(x1);

  return 0;
}
int average(int n, float *x, float *ave, float *rms){
  int    i;
  float  xsum,xsqrsum;

  if(n<0){
    fprintf(stderr,"Negative array length %d illegal\n",n);
    return 1;
  }
  if(n==0){*ave=*rms=0.0; return 0;}

  xsum=xsqrsum=0.0;
  for(i=0;i<n;i++){
    xsum    += x[i];
    xsqrsum += x[i]*x[i];                
  }
  *ave  = xsum/(float)n;
  if(n>1)
    *rms  = sqrt((xsqrsum -n*(*ave)*(*ave))/(float)(n-1));
  else
    *rms  = 0.0;

  return 0;
}
