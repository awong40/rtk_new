/* main function file  ----------------------------------------------------
*-----------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
int i,j,ret;
ret=estpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
    
    if (!ret) fprintf(stderr,"%40s\r","");
    return ret;
}
