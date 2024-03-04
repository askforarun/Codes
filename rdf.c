#include "stdheader.h"
#include "getljcoeff.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"
#include "particle.h"

// gcc rdf.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubecentres.c particle.c delU.c tailcorr.c -lm -o rdf
//./a.out dump_0.45.lammpstrj paircoeff.txt 

//function declarations-----------
int getnatoms(FILE *fp); 
double energy_diff(Box box, Atom *atom, double cutoff, Particle *particle);
int gettimestep(FILE *fp);
double tailcorr (Atomparam *atomparam, Particle *particle, double rc, double vol);
//------------------------------------------------------------------------------


FILE *fp1,*fp2, *fp3,*fp4,*fp5;
	
        double volavg=0;
        int natoms,timestep,Nx=20,Ny=20,Nz=20,counter=0,counter_frame=0,i=0,j=0;
        double rc=9,delUtemp=0,exp_minusbetadelUtemp=0,exp_minusbetadelUframe=0,exp_minusbetadelU=0;
        double delUtail,dx=0,dy=0,dz=0;
        double vol,void_vol;Particle *water;
        Box box; Atom *atom; Centres *centre; Atomparam *atomparam; Cell *cell;
        Atom *atom1; Atom *atom2;
       	char str1[100], str2[100];
	size_t nbytes = 100;
	char *my_string=NULL;
        double r;
        int *ind1, *ind2,*ind3,*ind4,totindex1=0,totindex2=0,ret=0;
        size_t len=0;
        ssize_t read;
        int getnum;
        int k=0;

 int bin=0;
double delr=0.1;
double dist=15;
int N;
double *H;



void main(int argc, char *argv[]) {
        
        fp1 = fopen(argv[1], "r");
        fp2=fopen(argv[2], "r");
        fp3=fopen(argv[3], "r");
        fp4=fopen(argv[4], "r");	


 N= (int) (dist/delr);
H=malloc(sizeof(double)*(N+1));
for(i=0;i<=N ;i++){
H[i]=0;
}


//printf("%d\n",N);
natoms=getnatoms(fp1);//Get the number of atoms

                                atom1=(Atom *) malloc(sizeof(Atom));
                                atom1->id = (int *) malloc(sizeof(int) *natoms);
                                atom1->type =(int *) malloc(sizeof(int) * natoms);
				atom1->x = (double *) malloc(sizeof(double) * natoms);
				atom1->y = (double *) malloc(sizeof(double) * natoms);
				atom1->z = (double *) malloc(sizeof(double) * natoms);
				atom1->q = (double *) malloc(sizeof(double) * natoms);
                                atom1->eps=(double *) malloc(sizeof(double) * natoms);
                                atom1->sigma=(double *) malloc(sizeof(double) * natoms);

                                 atom2=(Atom *) malloc(sizeof(Atom));
                                atom2->id = (int *) malloc(sizeof(int) *natoms);
                                atom2->type =(int *) malloc(sizeof(int) * natoms);
				atom2->x = (double *) malloc(sizeof(double) * natoms);
				atom2->y = (double *) malloc(sizeof(double) * natoms);
				atom2->z = (double *) malloc(sizeof(double) * natoms);
				atom2->q = (double *) malloc(sizeof(double) * natoms);
                                atom2->eps=(double *) malloc(sizeof(double) * natoms);
                                atom2->sigma=(double *) malloc(sizeof(double) * natoms);
                               
ind1=malloc(sizeof(int)*natoms); 
ind2=malloc(sizeof(int)*natoms);   
ind3=malloc(sizeof(int)*natoms);  
ind4=malloc(sizeof(int)*natoms);  

i=0; 
while ((read = getline(&my_string, &nbytes, fp3)) != -1) {
sscanf(my_string,"%d %d\n",&ind1[i],&ind2[i]);
//printf("%d %d %d\n",i,ind1[i],ind2[i]);
i=i+1;
}
totindex1=i;



i=0;
while ((read = getline(&my_string, &nbytes, fp4)) != -1) {
sscanf(my_string,"%d %d\n",&ind3[i],&ind4[i]);
//printf("%d %d %d\n",i,ind3[i],ind4[i]);
i=i+1;
}
totindex2=i;
     



//printf("%d %d\n",totindex1, totindex2);


//printf("%d %d\n",totindex1,totindex2);
atomparam=getparameters(fp1,fp2,natoms);

Ljcoeff *lj=Ljcoeffs(fp1,atomparam,natoms);//Get the Ljcoefficients epsilon and sigma for each atom in the system
//water=particleinfo();// Get epsilon, sigma, and coordinates,types of the atoms of the inserted particle (at equilibrium)-------------

rewind(fp1);//Start from begining of the dump file
//read each line of the dump file
//while (!feof(fp1)) {
while((read=getline(&my_string, &nbytes, fp1))!=-1) { 
		
		//getline(&my_string, &nbytes, fp1);
		//sscanf(my_string, "%s %s", str1, str2);
if(sscanf(my_string, "%s %s", str1, str2)==2){

//------------------------------------------------------------------------------               
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "TIMESTEP") == 0) {
        counter=counter+1;
	timestep=gettimestep(fp1);
     //  printf("%d\n",timestep);
}
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "BOX") == 0) {
              i=i+1;
              char boxfile[60];
             
             
              box=boxdim(fp1);
             
              vol= (box.xhi-box.xlo)*(box.yhi-box.ylo)*(box.zhi-box.zlo);
         //  printf("%lf\n",vol); 
              }
//----------------Widom insertion begins here ------------------------------------------------------

if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {
//if(timestep==0){

 printf("%d\n",timestep);

     // char filename1[60];  
 //sprintf(filename1, "box_%d.txt", counter);
// fp4 = fopen(filename1,"w");
 //fprintf(fp4,"%lf %lf %lf %lf %lf %lf", box.xlo,box.xhi,box.ylo,box.yhi,box.zlo,box.zhi);
//printf("%d\n",counter);
///delUtail=tailcorr(atomparam,water,rc,vol);
//printf("%lf\n",delUtail);
atom=atominfo(fp1,box,lj,natoms);


for(j=0;j<totindex1;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind1[j]){
atom1->x[j]=atom->x[i];
atom1->y[j]=atom->y[i];
atom1->z[j]=atom->z[i];
atom1->id[j]=atom->id[i];
//printf("%d %d %d\n",ind1[j],atom1->id[j],atom->id[i]);
}
//printf("%d %d\n",ind1[j],atom->id[i]);
}
}



for(j=0;j<totindex2;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind3[j]){
atom2->x[j]=atom->x[i];
atom2->y[j]=atom->y[i];
atom2->z[j]=atom->z[i];
atom2->id[j]=atom->id[i];
//printf("%d %d %d\n",ind3[j],atom2->id[j],atom->id[i]);
}
//printf("%d %d\n",ind1[j],atom->id[i]);
}
}

for (i = 0; i < totindex1; i++) {
    for (j = 0; j < totindex2; j++) {
        if(ind2[i]!=ind4[j] && i!=j){
                                     
                                        dx=atom1->x[i]-atom2->x[j];
                                        dy=atom1->y[i]-atom2->y[j];
                                        dz=atom1->z[i]-atom2->z[j];
                                        if (dx <= 0.5*(box.xhi-box.xlo)) {dx = dx+(box.xhi-box.xlo);}
                                        if (dx > 0.5*(box.xhi-box.xlo)) {dx = dx-(box.xhi-box.xlo);}
                                        if (dy <= 0.5*(box.yhi-box.ylo)) {dy = dy+(box.yhi-box.ylo);}
                                        if (dy > 0.5*(box.yhi-box.ylo)) {dy = dy-(box.yhi-box.ylo);}
                                        if (dz <= 0.5*(box.zhi-box.zlo)) {dz = dz+(box.zhi-box.zlo);}
                                        if (dz > 0.5*(box.zhi-box.zlo)) {dz = dz-(box.zhi-box.zlo);}
                                        r=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
                                         if(r < dist){
                                       
                                        bin=ceil(r/delr);
                                        H[bin]=H[bin]+1;
                                       //  printf("%d %lf %d %d %d %d\n",ind2[i],H[bin],i,j,atom1->id[i],atom2->id[j]);
                                        //printf("%d\n",bin);
                                       }
                                   



}
}
}





//autocorrelation code begins here. 

volavg=volavg+vol;
//printf("%lf %lf\n",volavg,vol);


free_atom(atom);
exp_minusbetadelU=exp_minusbetadelU+exp_minusbetadelUframe;
exp_minusbetadelUframe=0;
counter_frame=0;
}
//}
//---------------------Widom insertion ends here---------------------

   }




            
	}

char filename[60];  
sprintf(filename, "file_rdf.txt");
fp5 = fopen(filename,"w"); 
volavg=volavg/counter;
double density=totindex2/volavg;
for(i=1;i<=N ;i++){
r=i*delr;
H[i]=H[i]/(totindex1*counter*4*3.14*pow(r,2)*delr);
H[i]=H[i]/density;
fprintf(fp5,"%lf %lf\n",r,H[i]);
//printf("%lf %lf\n",r,H[i]);
}
         // printf("%lf %d\n",-log(exp_minusbetadelU/counter)+2*ecor,counter);
        // free_getljcoeff(lj); 
         //free(my_string);
         //free_getparameters(atomparam);
         fclose(fp1);
         fclose(fp2);
//         fclose(fp3);
//         flcose(fp4);
       
  }
