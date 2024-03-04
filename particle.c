#include "stdheader.h"
#include "particle.h"

Particle *particleinfo()
{
Particle *particle;
particle=(Particle *)malloc(sizeof(Particle));
particle->x=(double *)malloc(sizeof(double)*particle->natoms);
particle->y=(double *)malloc(sizeof(double)*particle->natoms);
particle->z=(double *)malloc(sizeof(double)*particle->natoms);
particle->natoms=1;
particle->lj=(LJcoeffparticle *)malloc(sizeof(LJcoeffparticle));
particle->lj->eps=(double *)malloc(sizeof(double)*particle->natoms);
particle->lj->sigma=(double *)malloc(sizeof(double)*particle->natoms);
particle->lj->eps[0]=1;
particle->lj->sigma[0]=1;
particle->x[0]=1;
particle->y[0]=1;
particle->z[0]=1;
return particle;
}



Particle *createparticle(Particle *particle, double x, double y, double z)
{

particle->x[0]=x;
particle->y[0]=y;
particle->z[0]=z;
//printf("%lf\n",particle->x[0]);
return particle;
}


void free_particle(Particle *particle)
{
free(particle->x);
free(particle->y);
free(particle->z);
free(particle->lj->eps);
free(particle->lj->sigma);
free(particle->lj);
free(particle);
}
