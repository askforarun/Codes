#include "stdheader.h"

size_t nbytes;
char *my_string=NULL;
int i,j,k;
double **dx,**dy,**dz;
double P2temp,**P2,*P2avg;
FILE *fp1;
ssize_t read;
char filename[100];
int Tlags;
char *tlags;

int numlines;
char *t;
int T;
void main(int argc, char *argv[]){ 



t= argv[1];
tlags = argv[2];
T = atoi(t);
Tlags = atoi(tlags);

fp1 = fopen("P2_1.txt","r");   

if (fp1 == NULL) {
    printf("%s\n",filename);
     exit(1);
}





while ((read = getline(&my_string, &nbytes, fp1))!= -1) {
j=j+1;
}
numlines=j;
//printf("%d\n",numlines);

dx    = malloc(numlines* sizeof(double *)); 
dy    = malloc(numlines* sizeof(double *)); 
dz    = malloc(numlines* sizeof(double *)); 
P2    = malloc(numlines* sizeof(double *)); 
P2avg = malloc(numlines* sizeof(double)); 

for (i=0; i<T; i++) {
     dx[i]= malloc(T* sizeof(double)); 
     dy[i]= malloc(T* sizeof(double)); 
     dz[i]= malloc(T* sizeof(double)); 
     P2[i]= malloc(T* sizeof(double)); 
      }


for(i=1;i<=T;i++){
sprintf(filename,"P2_%d.txt",i);
fp1 = fopen(filename,"r");   

if (fp1 == NULL) {
    // printf("%s\n",filename);
     exit(1);
}

j=0;
while ((read = getline(&my_string, &nbytes, fp1))!= -1) {
sscanf(my_string,"%lf %lf %lf",&dx[j][i-1],&dy[j][i-1],&dz[j][i-1]);
//printf("%lf %lf %lf\n",dx[j][i-1],dy[j][i-1],dz[j][i-1]); 
j=j+1;
}

fclose(fp1);
}


//--------------autocorrelation code begins here-----------

for(j=0;j<numlines;j++){
  for(i=0;i < Tlags;i++) {
           for(k=0; k < T-i;k++)  
          {
       P2temp=P2temp+dx[j][k]*dx[j][k+i]+dy[j][k]*dy[j][k+i]+dz[j][k]*dz[j][k+i];
          }
P2[j][i]=P2temp/(k+1);
//printf("%lf\n", P2[j][i]);
P2temp=0;
     }
}



for(i=0;i< Tlags;i++){
   for(j=0;j<numlines;j++){
P2avg[i] =  P2avg[i] + P2[j][i];
                      }
P2avg[i] =  P2avg[i]/numlines;
printf("%lf\n",P2avg[i]/P2avg[0]); //divide by variance.
}
//-----------------------------------------------------------



}



