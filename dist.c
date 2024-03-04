#include "stdheader.h"
#include "getljcoeff.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"


// gcc dist.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c -lm -o dist
//./a.out dump_550.lammpstrj allcoeffs_equalreactivity.txt P1-P1-1.txt P1-P1-2.txt 

//function declarations-----------
int getnatoms(FILE *fp); 
int gettimestep(FILE *fp,size_t nbytes,char *my_string);
//------------------------------------------------------------------------------


FILE *fp1,*fp2, *fp3,*fp4,*fp5,*fp6;
	
        int natoms,timestep,counter=0,i=0,j=0;
     
        double dx=0,dy=0,dz=0;
        double vol,void_vol;
        Box box; Atom *atom; Centres *centre; Atomparam *atomparam; Cell *cell;
        Atom *atom1; Atom *atom2;
       	char str1[100], str2[100];
	size_t nbytes = 100;
	char *my_string=NULL;
        char *line=NULL;
        double r;
        int *ind1, *ind2,totindex1=0,totindex2=0,ret=0;
        size_t len=0;
        ssize_t read;
        int getnum;
        char filename1[60];  
        char filename2[60];
        char boxfile[60];

void main(int argc, char *argv[]) {
fp1 =fopen(argv[1], "r");
fp2 =fopen(argv[2], "r");
fp3 =fopen(argv[3], "r");

	
  
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

while ((read = getline(&my_string, &nbytes, fp3)) != -1) {
sscanf(my_string,"%d %d",&ind1[i],&ind2[i]);
//printf("%d %d %d\n",i,ind1[i],ind2[i]);
i=i+1;
}

totindex1=i;

printf("%d\n",totindex1);

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
	timestep=gettimestep(fp1,nbytes,my_string);
       //printf("%d\n",timestep);
}
if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "BOX") == 0) {
      
      
              box=boxdim(fp1);
            
              vol= (box.xhi-box.xlo)*(box.yhi-box.ylo)*(box.zhi-box.zlo);
         //  printf("%lf\n",vol); 
              }


if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {
if(timestep <= 10000000){
 

      sprintf(boxfile, "box_%d.txt", counter);
  //            fp4 = fopen(boxfile,"w"); 
  //            fprintf(fp4,"%lf %lf %lf %lf %lf %lf" , box.xlo,box.xhi,box.ylo,box.yhi,box.zlo,box.zhi);
//printf("%lf %lf %lf %lf %lf %lf" , box.xlo,box.xhi,box.ylo,box.yhi,box.zlo,box.zhi);
     
    // fclose(fp4);
                   
 sprintf(filename1, "P2_%d.txt", counter);
 fp5 = fopen(filename1,"w");   
          
// sprintf(filename2, "file_%d", counter);
// fp6 = fopen(filename2,"w");   
atom=atominfo(fp1,box,lj,natoms);
//for(j=0;j<atom->natoms;j++){
//fprintf(fp6, "%d %lf %lf %lf\n", atom->id[j],atom->x[j],
//							atom->y[j], atom->z[j]);
//}

printf("%d\n",timestep);

///delUtail=tailcorr(atomparam,water,rc,vol);
//printf("%lf\n",delUtail);

   // fclose(fp6);

for(j=0;j<totindex1;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind1[j]){
atom1->x[j]=atom->x[i];
atom1->y[j]=atom->y[i];
atom1->z[j]=atom->z[i];
atom1->id[j]=atom->id[i];
//printf("%d %lf %lf\n",j,atom1->z[j],atom->z[i]);
}
}
}


//for(j=0;j<totindex1;j++){
for(j=0;j<totindex1;j++){
for (i = 0; i < atom->natoms; i++){
if(atom->id[i]==ind2[j]){
atom2->x[j]=atom->x[i];
atom2->y[j]=atom->y[i];
atom2->z[j]=atom->z[i];
atom2->id[j]=atom->id[i];
//printf("%d %lf %lf\n",j,atom2->z[j],atom->z[i]);
}
}
}




for (i = 0; i < totindex1; i++) {
           
           if (atom1->x[i]==0 || atom2->x[i]==0 || atom1->y[i]==0 || atom2->y[i]==0 || atom1->z[i]==0 || atom2->z[i]==0) {printf("%d %lf %lf %lf\n",i,atom1->x[i],atom1->y[i],atom1->z[i]);printf("%d %lf %lf %lf\n",i,atom2->x[i],atom2->y[i],atom2->z[i]);exit(1);}
                                        dx=atom1->x[i]-atom2->x[i];
                                        dy=atom1->y[i]-atom2->y[i];
                                        dz=atom1->z[i]-atom2->z[i];
                                        if (dx <= -0.5*(box.xhi-box.xlo)) {dx = dx+(box.xhi-box.xlo);}
                                        if (dx >   0.5*(box.xhi-box.xlo)) {dx = dx-(box.xhi-box.xlo);}
                                        if (dy <= -0.5*(box.yhi-box.ylo)) {dy = dy+(box.yhi-box.ylo);}
                                        if (dy >   0.5*(box.yhi-box.ylo)) {dy = dy-(box.yhi-box.ylo);}
                                        if (dz <= -0.5*(box.zhi-box.zlo)) {dz = dz+(box.zhi-box.zlo);}
                                        if (dz >   0.5*(box.zhi-box.zlo)) {dz = dz-(box.zhi-box.zlo);}

                                        r=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
                                      //  if (r >35) {
                                          //printf("%d %d %d %lf\n",i,atom1->id[i],atom2->id[i],r);exit(1);
                                        //  printf("%d %lf %lf %lf\n",atom1->id[i],atom1->x[i],atom1->y[i],atom1->z[i]);
                                         // printf("%d %lf %lf %lf\n",atom2->id[i],atom2->x[i],atom2->y[i],atom2->z[i]);exit(1);
                                         // }
                                     // printf("%lf\n",r);
                                      //printf(" %d %lf %lf %lf\n",i,dx,dy,dz);
                                        fprintf(fp5,"%lf %lf %lf %lf\n",dx,dy,dz,r);
                                       // fprintf(fp6,"%lf\n",r);
                                            
                               
		        	}



                                        fclose(fp5);    



free_atom(atom);

}
}


   }


            
	}
         
         fclose(fp1);
         fclose(fp2);

       
  }
