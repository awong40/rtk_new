/* main function file  ----------------------------------------------------
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"
int main(int argc, char **argv)
{

    prcopt_t opt_=opt_;
    double *rs,*dts,*var,*azel_,*resp;
    int i,n,stat,vsat,svh;




     rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
     if (n<=0) {

        return -2;
    }
    stat=estpos(n,rs,dts,var,svh,opt_,azel_,vsat,resp);
    if (!stat) fprintf(stderr,"%40s\r","");
    return stat;
}   
     
   
    
  
