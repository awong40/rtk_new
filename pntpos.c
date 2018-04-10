//
// Created by Anissa Wong on 4/10/18.
//

//
// Created by Anissa Wong on 4/10/18.
//

#include "rtklib.h"
#include <math.h>
static const char rcsid[]="$Id:$";

/* constants -----------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#define NX          (4+3)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */
const prcopt_t prcopt_default={ /* defaults processing options */
        PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
        15.0*D2R,{{0,0}},           /* elmin,snrmask */
        0,1,1,1,                    /* sateph,modear,glomodear,bdsmodear */
        5,0,10,                     /* glomodear,maxout,minlock,minfix */
        0,0,0,0,                    /* estion,esttrop,dynamics,tidecorr */
        1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
        0,0,                        /* rovpos,refpos */
        {100.0,100.0},              /* eratio[] */
        {100.0,0.003,0.003,0.0,1.0}, /* err[] */
        {30.0,0.03,0.3},            /* std[] */
        {1E-4,1E-3,1E-4,1E-1,1E-2}, /* prn[] */
        5E-12,                      /* sclkstab */
        {3.0,0.9999,0.20},          /* thresar */
        0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
        30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
        {0},{0},{0},                /* baseline,ru,rb */
        {"",""},                    /* anttype */
        {{0}},{{0}},{0}             /* antdel,pcv,exsats */
};


/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}

static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P,lam_L1;
    int i,j,nv=0,sys,mask[4]={0};

    trace(3,"resprng : n=%d\n",n);

    for (i=0;i<3;i++) rr[i]=x[i]; dtr=x[3];

    ecef2pos(rr,pos);

    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;

        if (!(sys=satsys(obs[i].sat,NULL))) continue;

        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
            trace(2,"duplicated observation data %s sat=%2d\n",
                  time_str(obs[i].time,3),obs[i].sat);
            i++;
            continue;
        }

        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);

        /* design matrix */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);

        /* time system and receiver bias offset correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else mask[0]=1;

        vsat[i]=1; resp[i]=v[nv]; (*ns)++;

        /* error variance */
        var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;

        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<4;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}



/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;

    trace(3,"estpos  : n=%d\n",n);

    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);

    for (i=0;i<3;i++) x[i]=sol->rr[i];

    for (i=0;i<MAXITR;i++) {

        /* pseudorange residuals */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);

        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* weight by variance */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) x[j]+=dx[j];

        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* gal-gps time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* bds-gps time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(unsigned char)ns;
            sol->age=sol->ratio=0.0;



            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);

    free(v); free(H); free(var);

    return 0;
}