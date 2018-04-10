/*------------------------------------------------------------------------------
* rtkcmn.c : rtklib common functions
*
*          Copyright (C) 2007-2018 by T.TAKASU, All rights reserved.
*
* options : -DLAPACK   use LAPACK/BLAS
*           -DMKL      use Intel MKL
*           -DTRACE    enable debug trace
*           -DWIN32    use WIN32 API
*           -DNOCALLOC no use calloc for zero matrix
*           -DIERS_MODEL use GMF instead of NMF
*           -DDLL      built for shared library
*
* references :
*     [1] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [2] RTCA/DO-229C, Minimum operational performanc standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [3] M.Rothacher, R.Schmid, ANTEX: The Antenna Exchange Format Version 1.4,
*         15 September, 2010
*     [4] A.Gelb ed., Applied Optimal Estimation, The M.I.T Press, 1974
*     [5] A.E.Niell, Global mapping functions for the atmosphere delay at radio
*         wavelengths, Jounal of geophysical research, 1996
*     [6] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.00, November 28, 2007
*     [7] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [8] China Satellite Navigation Office, BeiDou navigation satellite system
*         signal in space interface control document, open service signal B1I
*         (version 1.0), Dec 2012
*     [9] J.Boehm, A.Niell, P.Tregoning and H.Shuh, Global Mapping Function
*         (GMF): A new empirical mapping function base on numerical weather
*         model data, Geophysical Research Letters, 33, L07304, 2006
*     [10] GLONASS/GPS/Galileo/Compass/SBAS NV08C receiver series BINR interface
*         protocol specification ver.1.3, August, 2012
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/12 1.0 new
*           2007/03/06 1.1 input initial rover pos of pntpos()
*                          update only effective states of filter()
*                          fix bug of atan2() domain error
*           2007/04/11 1.2 add function antmodel()
*                          add gdop mask for pntpos()
*                          change constant MAXDTOE value
*           2007/05/25 1.3 add function execcmd(),expandpath()
*           2008/06/21 1.4 add funciton sortobs(),uniqeph(),screent()
*                          replace geodist() by sagnac correction way
*           2008/10/29 1.5 fix bug of ionosphereic mapping function
*                          fix bug of seasonal variation term of tropmapf
*           2008/12/27 1.6 add function tickget(), sleepms(), tracenav(),
*                          xyz2enu(), satposv(), pntvel(), covecef()
*           2009/03/12 1.7 fix bug on error-stop when localtime() returns NULL
*           2009/03/13 1.8 fix bug on time adjustment for summer time
*           2009/04/10 1.9 add function adjgpsweek(),getbits(),getbitu()
*                          add function geph2pos()
*           2009/06/08 1.10 add function seph2pos()
*           2009/11/28 1.11 change function pntpos()
*                           add function tracegnav(),tracepeph()
*           2009/12/22 1.12 change default parameter of ionos std
*                           valid under second for timeget()
*           2010/07/28 1.13 fix bug in tropmapf()
*                           added api:
*                               obs2code(),code2obs(),cross3(),normv3(),
*                               gst2time(),time2gst(),time_str(),timeset(),
*                               deg2dms(),dms2deg(),searchpcv(),antmodel_s(),
*                               tracehnav(),tracepclk(),reppath(),reppaths(),
*                               createdir()
*                           changed api:
*                               readpcv(),
*                           deleted api:
*                               uniqeph()
*           2010/08/20 1.14 omit to include mkl header files
*                           fix bug on chi-sqr(n) table
*           2010/12/11 1.15 added api:
*                               freeobs(),freenav(),ionppp()
*           2011/05/28 1.16 fix bug on half-hour offset by time2epoch()
*                           added api:
*                               uniqnav()
*           2012/06/09 1.17 add a leap second after 2012-6-30
*           2012/07/15 1.18 add api setbits(),setbitu(),utc2gmst()
*                           fix bug on interpolation of antenna pcv
*                           fix bug on str2num() for string with over 256 char
*                           add api readblq(),satexclude(),setcodepri(),
*                           getcodepri()
*                           change api obs2code(),code2obs(),antmodel()
*           2012/12/25 1.19 fix bug on satwavelen(),code2obs(),obs2code()
*                           add api testsnr()
*           2013/01/04 1.20 add api gpst2bdt(),bdt2gpst(),bdt2time(),time2bdt()
*                           readblq(),readerp(),geterp(),crc16()
*                           change api eci2ecef(),sunmoonpos()
*           2013/03/26 1.21 tickget() uses clock_gettime() for linux
*           2013/05/08 1.22 fix bug on nutation coefficients for ast_args()
*           2013/06/02 1.23 add #ifdef for undefined CLOCK_MONOTONIC_RAW
*           2013/09/01 1.24 fix bug on interpolation of satellite antenna pcv
*           2013/09/06 1.25 fix bug on extrapolation of erp
*           2014/04/27 1.26 add SYS_LEO for satellite system
*                           add BDS L1 code for RINEX 3.02 and RTCM 3.2
*                           support BDS L1 in satwavelen()
*           2014/05/29 1.27 fix bug on obs2code() to search obs code table
*           2014/08/26 1.28 fix problem on output of uncompress() for tar file
*                           add function to swap trace file with keywords
*           2014/10/21 1.29 strtok() -> strtok_r() in expath() for thread-safe
*                           add bdsmodear in procopt_default
*           2015/03/19 1.30 fix bug on interpolation of erp values in geterp()
*                           add leap second insertion before 2015/07/01 00:00
*                           add api read_leaps()
*           2017/01/03 1.31 add leap second before 2017/1/1 00:00:00
*           2018/01/29 1.32 chanage api crc16() -> rtk_crc16()
*                           chanage api crc32() -> rtk_crc32()
*                           chanage api crc24q() -> rtk_crc24q()
*-----------------------------------------------------------------------------*/
#define _POSIX_C_SOURCE 199309
#include <stdarg.h>
#include <ctype.h>
#ifndef WIN32
#include <dirent.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "rtklib.h"

static const char rcsid[]="$Id: rtkcmn.c,v 1.1 2008/07/17 21:48:06 ttaka Exp ttaka $";

/* constants -----------------------------------------------------------------*/
/* fatal error 339---------------------------------------------------------------*/
#ifdef MKL
#define LAPACK
#define dgemm_      dgemm
#define dgetrf_     dgetrf
#define dgetri_     dgetri
#define dgetrs_     dgetrs
#endif
#ifdef LAPACK
extern void dgemm_(char *, char *, int *, int *, int *, double *, double *,
                   int *, double *, int *, double *, double *, int *);
extern void dgetrf_(int *, int *, double *, int *, int *, int *);
extern void dgetri_(int *, double *, int *, int *, double *, int *, int *);
extern void dgetrs_(char *, int *, int *, double *, int *, int *, double *,
                    int *, int *);
#endif

#ifdef IERS_MODEL
extern int gmf_(double *mjd, double *lat, double *lon, double *hgt, double *zd,
                double *gmfh, double *gmfw);
#endif
static void fatalerr(const char *format, ...)
{
    va_list ap;
    va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
    exit(-9);
}
/* satellite number to satellite system 380----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/

extern int satsys(int sat, int *prn)
{
    int sys=0;
    if (sat<=0||MAXSAT<sat) sat=0;
    else if (sat<=NSATGPS) {
        sys=SYS_GPS; sat+=MINPRNGPS-1;
    }
    else if ((sat-=NSATGPS)<=NSATGLO) {
        sys=SYS_GLO; sat+=MINPRNGLO-1;
    }
    else if ((sat-=NSATGLO)<=NSATGAL) {
        sys=SYS_GAL; sat+=MINPRNGAL-1;
    }
    else if ((sat-=NSATQZS)<=NSATCMP) {
        sys=SYS_CMP; sat+=MINPRNCMP-1;
    }
    else if ((sat-=NSATLEO)<=NSATSBS) {
        sys=SYS_SBS; sat+=MINPRNSBS-1;
    }
    else sat=0;
    if (prn) *prn=sat;
    return sys;
}
/* new matrix 733------------------------------------------------------------------
* allocate memory of matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
    double *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
        fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* new integer matrix 748----------------------------------------------------------
* allocate memory of integer matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
    int *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
        fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* inner product 795---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* euclid norm 808-----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}
/* multiply matrix (wrapper of blas dgemm) 870-------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C)
* args   : char   *tr       I  transpose flags ("N":normal,"T":transpose)
*          int    n,k,m     I  size of (transposed) matrix A,B
*          double alpha     I  alpha
*          double *A,*B     I  (transposed) matrix A (n x m), B (m x k)
*          double beta      I  beta
*          double *C        IO matrix C (n x k)
* return : none
*-----------------------------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;
    
    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
           &ldb,&beta,C,&n);
}
/* inverse of matrix 884-----------------------------------------------------------
* inverse of matrix (A=A^-1)
* args   : double *A        IO  matrix (n x n)
*          int    n         I   size of matrix A
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double *work;
    int info,lwork=n*16,*ipiv=imat(n,1);
    
    work=mat(lwork,1);
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if (!info) dgetri_(&n,A,&n,ipiv,work,&lwork,&info);
    free(ipiv); free(work);
    return info;
}
/* least square estimation 1020-----------------------------------------------------
* least square estimation by solving normal equation (x=(A*A')^-1*A*y)
* args   : double *A        I   transpose of (weighted) design matrix (n x m)
*          double *y        I   (weighted) measurements (m x 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *x        O   estmated parameters (n x 1)
*          double *Q        O   esimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q)
{
    double *Ay;
    int info;
    
    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}
/* time to calendar day/time 1216---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
    const int mday[]={ /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t.time/86400);
    sec=(int)(t.time-(time_t)days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
/* add time 1330--------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
extern gtime_t timeadd(gtime_t t, double sec)
{
    double tt;
    
    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}
/* time to string 1504--------------------------------------------------------------
* convert gtime_t struct to string
* args   : gtime_t t        I   gtime_t struct
*          char   *s        O   string ("yyyy/mm/dd hh:mm:ss.ssss")
*          int    n         I   number of decimals
* return : none
*-----------------------------------------------------------------------------*/
extern void time2str(gtime_t t, char *s, int n)
{
    double ep[6];
    
    if (n<0) n=0; else if (n>12) n=12;
    if (1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* get time string 1528-------------------------------------------------------------
* get time string
* args   : gtime_t t        I   gtime_t struct
*          int    n         I   number of decimals
* return : time string
* notes  : not reentrant, do not use multiple in a function
*-----------------------------------------------------------------------------*/
extern char *time_str(gtime_t t, int n)
{
    static char buff[64];
    time2str(t,buff,n);
    return buff;
}
/* transform ecef to geodetic postion 1627------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* debug trace functions 2690-----------------------------------------------------*/
#ifdef TRACE

static FILE *fp_trace=NULL;     /* file pointer of trace */
static char file_trace[1024];   /* trace file */
static int level_trace=0;       /* level of trace */
static unsigned int tick_trace=0; /* tick time at traceopen (ms) */
static gtime_t time_trace={0};  /* time at traceopen */
static lock_t lock_trace;       /* lock for trace */

static void traceswap(void)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];
    
    lock(&lock_trace);
    
    if ((int)(time2gpst(time      ,NULL)/INT_SWAP_TRAC)==
        (int)(time2gpst(time_trace,NULL)/INT_SWAP_TRAC)) {
        unlock(&lock_trace);
        return;
    }
    time_trace=time;
    
    if (!reppath(file_trace,path,time,"","")) {
        unlock(&lock_trace);
        return;
    }
    if (fp_trace) fclose(fp_trace);
    
    if (!(fp_trace=fopen(path,"w"))) {
        fp_trace=stderr;
    }
    unlock(&lock_trace);
}
extern void traceopen(const char *file)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];
    
    reppath(file,path,time,"","");
    if (!*path||!(fp_trace=fopen(path,"w"))) fp_trace=stderr;
    strcpy(file_trace,file);
    tick_trace=tickget();
    time_trace=time;
    initlock(&lock_trace);
}
extern void traceclose(void)
{
    if (fp_trace&&fp_trace!=stderr) fclose(fp_trace);
    fp_trace=NULL;
    file_trace[0]='\0';
}
extern void tracelevel(int level)
{
    level_trace=level;
}
extern void trace(int level, const char *format, ...)
{
    va_list ap;
    
    /* print error message to stderr */
    if (level<=1) {
        va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
    }
    if (!fp_trace||level>level_trace) return;
    traceswap();
    fprintf(fp_trace,"%d ",level);
    va_start(ap,format); vfprintf(fp_trace,format,ap); va_end(ap);
    fflush(fp_trace);
}
extern void tracet(int level, const char *format, ...)
{
    va_list ap;
    
    if (!fp_trace||level>level_trace) return;
    traceswap();
    fprintf(fp_trace,"%d %9.3f: ",level,(tickget()-tick_trace)/1000.0);
    va_start(ap,format); vfprintf(fp_trace,format,ap); va_end(ap);
    fflush(fp_trace);
}
extern void tracemat(int level, const double *A, int n, int m, int p, int q)
{
    if (!fp_trace||level>level_trace) return;
    matfprint(A,n,m,p,q,fp_trace); fflush(fp_trace);
}
extern void traceobs(int level, const obsd_t *obs, int n)
{
    char str[64],id[16];
    int i;
    
    if (!fp_trace||level>level_trace) return;
    for (i=0;i<n;i++) {
        time2str(obs[i].time,str,3);
        satno2id(obs[i].sat,id);
        fprintf(fp_trace," (%2d) %s %-3s rcv%d %13.3f %13.3f %13.3f %13.3f %d %d %d %d %3.1f %3.1f\n",
              i+1,str,id,obs[i].rcv,obs[i].L[0],obs[i].L[1],obs[i].P[0],
              obs[i].P[1],obs[i].LLI[0],obs[i].LLI[1],obs[i].code[0],
              obs[i].code[1],obs[i].SNR[0]*0.25,obs[i].SNR[1]*0.25);
    }
    fflush(fp_trace);
}
extern void tracenav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!fp_trace||level>level_trace) return;
    for (i=0;i<nav->n;i++) {
        time2str(nav->eph[i].toe,s1,0);
        time2str(nav->eph[i].ttr,s2,0);
        satno2id(nav->eph[i].sat,id);
        fprintf(fp_trace,"(%3d) %-3s : %s %s %3d %3d %02x\n",i+1,
                id,s1,s2,nav->eph[i].iode,nav->eph[i].iodc,nav->eph[i].svh);
    }
    fprintf(fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gps[0],
            nav->ion_gps[1],nav->ion_gps[2],nav->ion_gps[3]);
    fprintf(fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gps[4],
            nav->ion_gps[5],nav->ion_gps[6],nav->ion_gps[7]);
    fprintf(fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gal[0],
            nav->ion_gal[1],nav->ion_gal[2],nav->ion_gal[3]);
}
extern void tracegnav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!fp_trace||level>level_trace) return;
    for (i=0;i<nav->ng;i++) {
        time2str(nav->geph[i].toe,s1,0);
        time2str(nav->geph[i].tof,s2,0);
        satno2id(nav->geph[i].sat,id);
        fprintf(fp_trace,"(%3d) %-3s : %s %s %2d %2d %8.3f\n",i+1,
                id,s1,s2,nav->geph[i].frq,nav->geph[i].svh,nav->geph[i].taun*1E6);
    }
}
extern void tracehnav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!fp_trace||level>level_trace) return;
    for (i=0;i<nav->ns;i++) {
        time2str(nav->seph[i].t0,s1,0);
        time2str(nav->seph[i].tof,s2,0);
        satno2id(nav->seph[i].sat,id);
        fprintf(fp_trace,"(%3d) %-3s : %s %s %2d %2d\n",i+1,
                id,s1,s2,nav->seph[i].svh,nav->seph[i].sva);
    }
}
extern void tracepeph(int level, const nav_t *nav)
{
    char s[64],id[16];
    int i,j;
    
    if (!fp_trace||level>level_trace) return;
    
    for (i=0;i<nav->ne;i++) {
        time2str(nav->peph[i].time,s,0);
        for (j=0;j<MAXSAT;j++) {
            satno2id(j+1,id);
            fprintf(fp_trace,"%-3s %d %-3s %13.3f %13.3f %13.3f %13.3f %6.3f %6.3f %6.3f %6.3f\n",
                    s,nav->peph[i].index,id,
                    nav->peph[i].pos[j][0],nav->peph[i].pos[j][1],
                    nav->peph[i].pos[j][2],nav->peph[i].pos[j][3]*1E9,
                    nav->peph[i].std[j][0],nav->peph[i].std[j][1],
                    nav->peph[i].std[j][2],nav->peph[i].std[j][3]*1E9);
        }
    }
}
extern void tracepclk(int level, const nav_t *nav)
{
    char s[64],id[16];
    int i,j;
    
    if (!fp_trace||level>level_trace) return;
    
    for (i=0;i<nav->nc;i++) {
        time2str(nav->pclk[i].time,s,0);
        for (j=0;j<MAXSAT;j++) {
            satno2id(j+1,id);
            fprintf(fp_trace,"%-3s %d %-3s %13.3f %6.3f\n",
                    s,nav->pclk[i].index,id,
                    nav->pclk[i].clk[j][0]*1E9,nav->pclk[i].std[j][0]*1E9);
        }
    }
}
extern void traceb(int level, const unsigned char *p, int n)
{
    int i;
    if (!fp_trace||level>level_trace) return;
    for (i=0;i<n;i++) fprintf(fp_trace,"%02X%s",*p++,i%8==7?" ":"");
    fprintf(fp_trace,"\n");
}
#else
extern void traceopen(const char *file) {}
extern void traceclose(void) {}
extern void tracelevel(int level) {}
extern void trace   (int level, const char *format, ...) {}
extern void tracet  (int level, const char *format, ...) {}
extern void tracemat(int level, const double *A, int n, int m, int p, int q) {}
extern void traceobs(int level, const obsd_t *obs, int n) {}
extern void tracenav(int level, const nav_t *nav) {}
extern void tracegnav(int level, const nav_t *nav) {}
extern void tracehnav(int level, const nav_t *nav) {}
extern void tracepeph(int level, const nav_t *nav) {}
extern void tracepclk(int level, const nav_t *nav) {}
extern void traceb  (int level, const unsigned char *p, int n) {}

#endif /* TRACE */