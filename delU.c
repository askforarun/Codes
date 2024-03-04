#include "stdheader.h"
#include "Atom.h"
#include "getboxdim.h"
#include "particle.h"


double energy_diff(Box box, Atom *atom, double cutoff, Particle *particle) {
int i,j;
double delU=0,delUtemp,delUtail=0,delUtailtemp=0,dx=0,dy,dz,r,eps_ij,sigma_ij,pi=3.14159;
//This part of the code calculates the total potential in the system.
//---------------------------------------------------------------------------------------------------------------
for (i = 0; i < atom->natoms-1; i++) {
                             //  printf("%d\n",i);
				for (j = i+1; j < atom->natoms; j++) {
                                        dx=atom->x[i]-atom->x[j];
                                        dy=atom->y[i]-atom->y[j];
                                        dz=atom->z[i]-atom->z[j];
                                        if (dx <= 0.5*(box.xhi-box.xlo)) {dx = dx+(box.xhi-box.xlo);}
                                        if (dx > 0.5*(box.xhi-box.xlo)) {dx = dx-(box.xhi-box.xlo);}
                                        if (dy <= 0.5*(box.yhi-box.ylo)) {dy = dy+(box.yhi-box.ylo);}
                                        if (dy > 0.5*(box.yhi-box.ylo)) {dy = dy-(box.yhi-box.ylo);}
                                        if (dz <= 0.5*(box.zhi-box.zlo)) {dz = dz+(box.zhi-box.zlo);}
                                        if (dz > 0.5*(box.zhi-box.zlo)) {dz = dz-(box.zhi-box.zlo);}
                                        r=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
                                       if(r< cutoff){
                                        // printf("%d %d\n",i,j);
                                        sigma_ij = 0.5*(atom->sigma[i] + atom->sigma[j]);
					eps_ij = sqrt(atom->eps[i] * atom->eps[j]);
                                       delUtemp=4*eps_ij*((pow((sigma_ij/r),12)-pow((sigma_ij/r),6)));
                                      // delUtailtemp= (8/3)*pi*(natoms/vol)*eps_ij*pow(sigma_ij,3)*((1/3)*pow((sigma_ij/cutoff),9)-pow((sigma_ij/cutoff),3));
                                      // delUtail=  delUtail+delUtailtemp;
                                       delU=delU+delUtemp;
                                                    }
                                          
				                                }
		        	}
return delU;
}
//---------------------------------------------------------------------------------------------------------------------
