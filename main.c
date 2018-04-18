/* main function file  ----------------------------------------------------
*-----------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
int i,j,ret;
ret=estpos(obs,n,rs,dts,nav,&opt_,sol,azel,vsat,resp,msg);
/* ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"",""); */
    
    if (!ret) fprintf(stderr,"%40s\r","");
    return ret;
}
