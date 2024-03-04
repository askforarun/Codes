#include "stdheader.h"
#include "particle.h"
#include "getljcoeff.h"

double tailcorr (Atomparam *atomparam, Particle *particle, double rc, double vol)
{
int i,j;
double sigma_ij, eps_ij,delUtail=0,delUtailtemp=0,pi=3.14159;
for(i=1;i <= atomparam->maxtype;i++)
{
for(j=1; j <=atomparam->maxtype;j++)
{
sigma_ij = 0.5*(atomparam->sigma[i] + atomparam->sigma[j]);
eps_ij = sqrt(atomparam->eps[i] * atomparam->eps[j]);
delUtailtemp= 8*pi/3*atomparam->typecount[i]*atomparam->typecount[j]*eps_ij*pow(sigma_ij,3)*(pow((sigma_ij/rc),9)/3-pow((sigma_ij/rc),3));
delUtail=delUtail+delUtailtemp;
//printf("%lf\n",4*eps_ij*((pow((sigma_ij/cutoff),12)-pow((sigma_ij/cutoff),6))));
//printf("%d %d %lf\n",atomparam->typecount[j],atomparam->typecount[i],rc);
}
}

delUtail=delUtail/vol;
return delUtail;
}


