#include "stdheader.h"
#include "getnatoms.h"
#include "getljcoeff.h"
#include "gettimestep.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"
#include "particle.h"


double energy_diff(Box box, Atom *atom, int natoms, double cutoff, Ljcoeff *lj, Particle *particle, LJcoeffparticle *ljcoeffparticle);

int main(int argc, char *argv[]) {
	FILE *fp1,*fp2;
	fp1 = fopen(argv[1], "r");
        fp2=fopen(argv[2], "r");
	int natoms,timestep,Nx=20,Ny=20,Nz=20,i,counter=0,counter_frame=0;
        double rc=2,delUtemp=0,exp_minusbetadelUtemp=0,exp_minusbetadelUframe=0,exp_minusbetadelU=0,epsilon=0;
        double T=1.2,void_vol,ecor;
        Box box; Atom *atom;
        Particle *particle; LJcoeffparticle *ljcoeffparticle;
       	char str1[100], str2[100];
	size_t nbytes = 100;
	char *my_string;
//-------------epsilon and sigma of the inserted particle------------
ljcoeffparticle=malloc(sizeof(LJcoeffparticle));
ljcoeffparticle->eps=malloc(sizeof(double));
ljcoeffparticle->sigma=malloc(sizeof(double));
ljcoeffparticle->eps[0]=1;
ljcoeffparticle->sigma[0]=1;
ljcoeffparticle->natoms=1;
//----------------------------

ecor=(8/3)*3.14159*0.45*((1/3)*(1/pow(rc,9))-(1/pow(rc,3)));

natoms=getnatoms(fp1);//Get the number of atoms
Ljcoeff *lj=Ljcoeffs(fp1,fp2,natoms);//Get the Ljcoefficients epsilon and sigma for each atom in the system
rewind(fp1);

//read each line of the dump file
while (!feof(fp1)) {
		my_string = (char *) malloc(nbytes + 1);
		getline(&my_string, &nbytes, fp1);
		sscanf(my_string, "%s %s", str1, str2);
if(sscanf(my_string, "%s %s", str1, str2)==2){
//------------------------------------------------------------------------------               
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "TIMESTEP") == 0) {
        counter=counter+1;
	timestep=gettimestep(fp1,nbytes,my_string,str1,str2);
       // printf("%d\n",timestep);
}
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "BOX") == 0) {
             box=boxdim(fp1,nbytes,my_string,str1,str2);
           //  printf("%lf %lf\n",(box.xhi-box.xlo)*(box.yhi-box.ylo)*(box.zhi-box.zlo),(double) Nx*Ny*Nz); 
              }
//-----------------------------------------------------------------------

if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {
atom=atominfo(fp1,nbytes,my_string,str1,str2,natoms);
Centres *centre=cubecentres(box,Nx,Ny,Nz);

for (i=0; i<(Nx*Ny*Nz);i++){
Particle *particle=createparticle(centre->xc[i], centre->yc[i], centre->zc[i]);
delUtemp=energy_diff(box, atom, natoms, rc, lj, particle, ljcoeffparticle);
exp_minusbetadelUtemp=exp(-(1/T)*delUtemp);
//printf("%d %lf %lf\n",i,exp_minusbetadelUtemp,delUtemp);
if(exp_minusbetadelUtemp >= epsilon){
counter_frame=counter_frame+1;

exp_minusbetadelUframe=exp_minusbetadelUframe+exp_minusbetadelUtemp;
                                    }
                           };

void_vol=(double)(counter_frame/(double)(Nx*Ny*Nz));
exp_minusbetadelUframe=exp_minusbetadelUframe/(Nx*Ny*Nz);
printf("%d %lf\n",counter,void_vol);

free_atom(atom);
free_cubecentres(centre);
}

exp_minusbetadelU=exp_minusbetadelU+exp_minusbetadelUframe;
exp_minusbetadelUframe=0;
counter_frame=0;

   }

//---------------------Widom insertion ends here---------------------
            
	}
        printf("%lf %d\n",-log(exp_minusbetadelU/counter)+2*ecor,counter);
        free_getljcoeff(lj);
        fclose(fp1);
	return 0;
}
