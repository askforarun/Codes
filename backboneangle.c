#include "stdheader.h" 
#include "getljcoeff.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"
#include "particle.h"

// gcc backboneangle.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubecentres.c particle.c delU.c tailcorr.c -lm -o bb
//./a.out $var3/dump_${k}.lammpstrj  $var4/allcoeffs_equalreactivity.txt index${i} temperature string

// index file contains 3 columns
//function declarations-----------
int getnatoms(FILE *fp); 
double energy_diff(Box box, Atom *atom, double cutoff, Particle *particle);
int gettimestep(FILE *fp);
double tailcorr (Atomparam *atomparam, Particle *particle, double rc, double vol);
//------------------------------------------------------------------------------

        FILE *fp1,*fp2, *fp3,*fp4;
        int fp5;
         
         int temp;
        double volavg=0;
        int natoms,timestep,Nx=20,Ny=20,Nz=20,counter=0,counter_frame=0,i=0,j=0;
        double rc=9,delUtemp=0,exp_minusbetadelUtemp=0,exp_minusbetadelUframe=0,exp_minusbetadelU=0;
        double delUtail,dx1=0,dy1=0,dz1=0,dx2=0,dy2=0,dz2=0,dx3=0,dy3=0,dz3=0;
        double vol,void_vol;Particle *water;
        double num=0,den=0,r1=0,r2=0,r3=0,theta=0;
        Box box; Atom *atom; Centres *centre; Atomparam *atomparam; Cell *cell;
        Atom *atom1, *atom2, *atom3;
       	char str1[100], str2[100];
	size_t nbytes = 100;
	char *my_string = NULL;
        char *line=NULL;
        double r;
        int *ind1, *ind2,*ind3,totindex1=0,totindex2=0,totindex3=0,ret=0,t;
        size_t len=0;
        ssize_t read;
        int getnum;
        int k=0;
 
char filename[60];  
 int bin=0;
double delr=0.1;
double dist=15;
int N;
double *H;


void main(int argc, char *argv[]) {

        fp1 = fopen(argv[1], "r");
        fp2=fopen(argv[2], "r");
        fp3=fopen(argv[3], "r");
        
        temp = atoi(argv[4]);
      
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



                                atom3=(Atom *) malloc(sizeof(Atom));
                                atom3->id = (int *) malloc(sizeof(int) *natoms);
                                atom3->type =(int *) malloc(sizeof(int) * natoms);
				atom3->x = (double *) malloc(sizeof(double) * natoms);
				atom3->y = (double *) malloc(sizeof(double) * natoms);
				atom3->z = (double *) malloc(sizeof(double) * natoms);
				atom3->q = (double *) malloc(sizeof(double) * natoms);
                                atom3->eps=(double *) malloc(sizeof(double) * natoms);
                                atom3->sigma=(double *) malloc(sizeof(double) * natoms);

                               
ind1=malloc(sizeof(int)*natoms);  
ind2=malloc(sizeof(int)*natoms);  
ind3=malloc(sizeof(int)*natoms);  
i=0; 

while ((read = getline(&my_string, &nbytes, fp3)) != -1) {
sscanf(my_string,"%d %d %d",&ind1[i],&ind2[i],&ind3[i]);
//printf("%d %d %d %d\n",i,ind1[i],ind2[i],ind3[i]);
i=i+1;
}
totindex1=i;




//printf("%d\n",totindex1);


//printf("%d %d\n",totindex1,totindex2);
atomparam=getparameters(fp1,fp2,natoms);

Ljcoeff *lj=Ljcoeffs(fp1,atomparam,natoms);//Get the Ljcoefficients epsilon and sigma for each atom in the system
//water=particleinfo();// Get epsilon, sigma, and coordinates,types of the atoms of the inserted particle (at equilibrium)-------------

rewind(fp1);//Start from begining of the dump file
//read each line of the dump file
while ((read = getline(&my_string, &nbytes, fp1)) != -1) {
		
		
		//sscanf(my_string, "%s %s", str1, str2);
if(sscanf(my_string, "%s %s", str1, str2)==2){

//------------------------------------------------------------------------------               
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "TIMESTEP") == 0) {
        counter=counter+1;
	timestep=gettimestep(fp1);
      
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

if (timestep >=0 && timestep <=1000000) {
printf("%d\n",timestep);
sprintf(filename, "dist_%d_%d_%s",counter,temp,argv[5]);
fp4=fopen(filename,"w");
//printf("%d\n",counter);



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


for(j=0;j<totindex1;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind2[j]){
atom2->x[j]=atom->x[i];
atom2->y[j]=atom->y[i];
atom2->z[j]=atom->z[i];
atom2->id[j]=atom->id[i];
//printf("%d %d %d\n",ind2[j],atom2->id[j],atom->id[i]);
}
//printf("%d %d\n",ind1[j],atom->id[i]);
}
}


for(j=0;j<totindex1;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind3[j]){
atom3->x[j]=atom->x[i];
atom3->y[j]=atom->y[i];
atom3->z[j]=atom->z[i];
atom3->id[j]=atom->id[i];
//printf("%d %d %d\n",ind3[j],atom3->id[j],atom->id[i]);
}
//printf("%d %d\n",ind1[j],atom->id[i]);
}
}




for (i = 0; i < totindex1; i++) {
 // for (j = i+1; j < totindex1; j++) {
//if(atom1->id[i]==107 && atom2->id[i]==156 && atom3->id[i]==130){
 if (atom1->x[i]==0 || atom2->x[i]==0 || atom3->x[i]==0 || atom1->y[i]==0 || atom2->y[i]==0 || atom3->y[i]==0 || atom1->z[i]==0 || atom2->z[i]==0 || atom3->z[i]==0) {printf("%d %s",i,"Atom ids not assigned properly\n");exit(1);}
                                        dx1=atom1->x[i]-atom2->x[i];
                                        dy1=atom1->y[i]-atom2->y[i];
                                        dz1=atom1->z[i]-atom2->z[i];
                                        dx2=atom1->x[i]-atom3->x[i];
                                        dy2=atom1->y[i]-atom3->y[i];
                                        dz2=atom1->z[i]-atom3->z[i];

                                        dx3=atom2->x[i]-atom3->x[i];
                                        dy3=atom2->y[i]-atom3->y[i];
                                        dz3=atom2->z[i]-atom3->z[i];
                                        if (dx3 <= -0.5*(box.xhi-box.xlo)) {dx3 = dx3+(box.xhi-box.xlo);}
                                        if (dx3 >   0.5*(box.xhi-box.xlo)) {dx3 = dx3-(box.xhi-box.xlo);}
                                        if (dy3 <= -0.5*(box.yhi-box.ylo)) {dy3 = dy3+(box.yhi-box.ylo);}
                                        if (dy3 >   0.5*(box.yhi-box.ylo)) {dy3 = dy3-(box.yhi-box.ylo);}
                                        if (dz3 <= -0.5*(box.zhi-box.zlo)) {dz3 = dz3+(box.zhi-box.zlo);}
                                        if (dz3 >   0.5*(box.zhi-box.zlo)) {dz3 = dz3-(box.zhi-box.zlo);}
                                        

                                        if (dx1 <= -0.5*(box.xhi-box.xlo)) {dx1 = dx1+(box.xhi-box.xlo);}
                                        if (dx1 >   0.5*(box.xhi-box.xlo)) {dx1 = dx1-(box.xhi-box.xlo);}
                                        if (dx2 <= -0.5*(box.xhi-box.xlo)) {dx2 = dx2+(box.xhi-box.xlo);}
                                        if (dx2 >   0.5*(box.xhi-box.xlo)) {dx2 = dx2-(box.xhi-box.xlo);}
                                        if (dy1 <= -0.5*(box.yhi-box.ylo)) {dy1 = dy1+(box.yhi-box.ylo);}
                                        if (dy1 >   0.5*(box.yhi-box.ylo)) {dy1 = dy1-(box.yhi-box.ylo);}
                                        if (dy2 <= -0.5*(box.yhi-box.ylo)) {dy2 = dy2+(box.yhi-box.ylo);}
                                        if (dy2 >   0.5*(box.yhi-box.ylo)) {dy2 = dy2-(box.yhi-box.ylo);}
                                        if (dz1 <= -0.5*(box.zhi-box.zlo)) {dz1 = dz1+(box.zhi-box.zlo);}
                                        if (dz1 >   0.5*(box.zhi-box.zlo)) {dz1 = dz1-(box.zhi-box.zlo);}
                                        if (dz2 <= -0.5*(box.zhi-box.zlo)) {dz2 = dz2+(box.zhi-box.zlo);}
                                        if (dz2 >   0.5*(box.zhi-box.zlo)) {dz2 = dz2-(box.zhi-box.zlo);}
                                        r1=sqrt(pow(dx1,2)+pow(dy1,2)+pow(dz1,2));
                                        r2=sqrt(pow(dx2,2)+pow(dy2,2)+pow(dz2,2));

                                        r3=sqrt(pow(dx3,2)+pow(dy3,2)+pow(dz3,2));
                                        num = dx1*dx2 + dy1*dy2 + dz1*dz2;
                                        den = r1*r2;
                                        theta = acos(num/den)*(180/3.1415);
                                        fprintf(fp4,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",i,atom1->id[i],atom2->id[i],atom3->id[i],dx1,dy1,dz1,dx2,dy2,dz2,r3,theta);
                                       // printf("%d %lf %lf %lf %lf %lf %lf %lf\n",i,dx1,dy1,dz1,dx2,dy2,dz2,theta);
                                      //autocorrelation code begins here. 

volavg=volavg+vol;




}
fclose(fp4);
}
}

//---------------------Widom insertion ends here---------------------

   }



            
	

   
}

       
  }
