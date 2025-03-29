# include <stdio.h>
# include <stdlib.h>
# include <errno.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <time.h>
# include "newcorr.h"
# include "svio.h"


void swap_bytes(unsigned short *p, int n)
{ while (n-- > 0) p[n] = p[n]/256 + (p[n] % 256)*256; }
void swap_short(unsigned *p, int n)
{ while (n-- > 0) p[n] = p[n]/65536 + (p[n] % 65536)*65536; }
void swap_long(void *p, int n)
{ swap_short(p, n); swap_bytes(p, 2*n); }
void swap_d(double *p, int n)
{ float temp, *v;
  swap_long(p, 2*n) ;
  while(n-- > 0)
  { v = (float *)p ;
    temp = *v; *v = v[1] ; v[1] = temp ;
    p++ ;
  }
} 
double lmst(double mjd)
/* gets local mean sidereal time for time_of_day given by mjd  */
{ double ut,tu,res,lmstime ;
  ut = mjd ;
  tu = (ut - 51544.5) / 36525.;  /* centuries since J2000  */
  res = ut + 74.05/360. ; res = res - floor(res) ;
  
  lmstime = res + (((0.093104 - tu * 6.2e-6) * tu
                  + 8640184.812866) * tu + 24110.54841)/86400.0 ;
  lmstime = (lmstime - floor(lmstime))*2.0 * M_PI ;
  return lmstime ;
}
char *mjd2iau_date( double mjd )
{
    double J = mjd + 2400000.5 ;  /*  is this correct ?  crs/june 2000 */
    static char date[256] ;
    int month, day;
    long year, a, c, d, x, y, jd;
    double dd;

    if( J < 1721425.5 ) return 0 ; /* January 1.0, 1 A.D.  dont accept BC ! */

    jd = J + 0.5; /* round Julian date up to integer */

    /* Find the number of Gregorian centuries
     * since March 1, 4801 B.C.
     */
    a = (100*jd + 3204500L)/3652425L;

    /* Transform to Julian calendar by adding in Gregorian century years
     * that are not leap years.
     * Subtract 97 days to shift origin of JD to March 1.
     * Add 122 days for magic arithmetic algorithm.
     * Add four years to ensure the first leap year is detected.
     */
    c = jd + 1486;
    if( jd >= 2299160.5 )
        c += a - a/4;
    else
        c += 38;
    /* Offset 122 days, which is where the magic arithmetic
     * month formula sequence starts (March 1 = 4 * 30.6 = 122.4).
     */
    d = (100*c - 12210L)/36525L;
    x = (36525L * d)/100L; /* Days in that many whole Julian years */

    /* Find month and day. */
    y = ((c-x)*100L)/3061L;
    day = c - x - ((306L*y)/10L);
    month = y - 1;
    if( y > 13 ) month -= 12;

    /* Get the year right. */
    year = d - 4715;
    if( month > 2 ) year -= 1;

    a = (jd + 1) % 7; /* Day of the week. */

    dd = day + J - jd + 0.5; /* Fractional part of day. */
    day = dd;
    { int h,m,s ;
      s = (dd-day)*24*3600 ;
      h = s/3600 ;
      m = (s-h*3600)/60 ;
      s -= (h*3600 + m*60) ;
      sprintf(date,"%04ld-%02d-%02dT%02d:%02d:%02d / UT",year,month,day,h,m,s) ;
    }
    return date ;
}


