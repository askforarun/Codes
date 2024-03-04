#include "stdheader.h"
#include "getljcoeff.h"

Atomparam *getparameters(FILE *fp1,FILE *fp2,int natoms)
{
int i=0,j=0,maxtype=0,count=0;
size_t nbytes = 100;
char *my_string =  malloc(nbytes);
char str1[100], str2[100];
Atomparam *atomparam;
atomparam=malloc(sizeof(Atomparam));
atomparam->eps = malloc(sizeof(double) * 1000);
atomparam->sigma = malloc(sizeof(double) *1000);

for (i = 1; i < 100; i++) {
getline(&my_string, &nbytes, fp2);
sscanf(my_string, "%s", str1);
if(sscanf(my_string, "%s", str1)==1){
if (strcmp(str1, "pair_coeff") == 0) {

//sscanf(my_string,"%*s %d",&j);
//sscanf(my_string,
//			"%*s %*d %*d %lf %lf",
//			&atomparam->eps[j], &atomparam->sigma[j]);



sscanf(my_string,
                        "%*s %*d %d %*[^/]%*[/]%*[^/]%*[/]%*[^/]%*[/]%*s",
                        &j);
sscanf(my_string,
                        "%*[^/]%*[/]%*[^/]%*[/]%*[^/]%*[/]%*s %lf %lf",
                        &atomparam->eps[j], &atomparam->sigma[j]);

//printf("%d %lf %lf\n",j,atomparam->eps[j],atomparam->sigma[j]);

}
}
}

atomparam->maxtype=j;
rewind(fp1);




atomparam->type= malloc(sizeof(int) *natoms);
rewind(fp1);

for (i = 0; i < 10; i++) {
getline(&my_string, &nbytes, fp1);
sscanf(my_string, "%s %s", str1, str2);
if(strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {
		for (j = 0; j < natoms; j++) {
                getline(&my_string, &nbytes, fp1);
		sscanf(my_string, "%*d %d", &atomparam->type[j]);
		//printf("%d\n",atomparam->type[j]);
             
	}
  break;
}
}

//#----------------------------
atomparam->typecount= malloc(sizeof(int) *natoms);
for(j=1;j<=atomparam->maxtype;j++){
for(i=0;i<natoms;i++){
if(atomparam->type[i]==j){
                count=count+1;
}
}

atomparam->typecount[j]=count;
//printf("%d %d\n",j,atomparam->typecount[j]);
count=0;
}
//--------------------------------------
//printf("%d\n",atomparam->maxtype);
free(my_string);
return atomparam; 
}


Ljcoeff *Ljcoeffs(FILE *fp1, Atomparam *atomparam, int natoms)
{
rewind(fp1);
int i=0,j=0;
size_t nbytes = 100;
char str1[100], str2[100];
double *eps_dummy, *sigma_dummy, *eps, *sigma;
char *my_string = malloc(nbytes);
Ljcoeff *ljcoeffs=malloc(sizeof(Ljcoeff));
ljcoeffs->eps =  malloc(sizeof(double) * natoms);
ljcoeffs->sigma =  malloc(sizeof(double) * natoms);
//loop of over all types, loop over all atoms and assign epsion and sigma to each atom based on the type of each atom
for(j=1;j<=atomparam->maxtype;j++){
for(i=0;i<natoms;i++){
if(atomparam->type[i]==j){
            
                ljcoeffs->eps[i] = atomparam->eps[j];
                ljcoeffs->sigma[i] = atomparam->sigma[j];
           
}
}
}
free(my_string);
return ljcoeffs;
}


void free_getparameters(Atomparam *atomparam)
{
free(atomparam->eps);
free(atomparam->sigma);
free(atomparam->type);
free(atomparam->typecount);
free(atomparam);
}


void free_getljcoeff(Ljcoeff *lj)
{
free(lj->eps);
free(lj->sigma);
free(lj);
}




