/* input ew matlab file NX */
/* column [0] = pseudorange residuals */
/* column [1] = obs bit stream */
/* column [2] = */
/* column [3] = */

#define SQR(x)      ((x)*(x))
#define NX          (4+3)       /* # of estimated parameters */
#define MAXITR      10          /* max number of iteration for point pos */


extern int estpos(const nv *nv, int n)
double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
int i,j,k,info,ns;
 v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,nv;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
/*    for (i=0;i<3;i++) x[i]=sol->rr[i]; */
    
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals */
/*        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp, */
/*                   &ns); */
        
        if (nv<NX) {
/*            sprintf(msg,"lack of valid sats ns=%d",nv); */
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
/*            sprintf(msg,"lsq error info=%d",info); */
            break;
        }

           
           

    }

}
