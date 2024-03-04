#include "stdheader.h"
#include "cubecentres.h"

Centres *cubecentres(Box box, int Nx,int Ny,int Nz){
  
                        int i,j,k,count=0;
                        double delx,dely,delz,x1,y1,z1,x2,y2,z2,dx,dy,dz;

                         Centres *centres;
                         centres=(Centres *)malloc(sizeof(Centres));
                         centres->xc=(double *)malloc(sizeof(double)*Nx*Ny*Nz);
                         centres->yc=(double *)malloc(sizeof(double)*Nx*Ny*Nz);
                         centres->zc=(double *)malloc(sizeof(double)*Nx*Ny*Nz);
                         delx=(box.xhi-box.xlo)/Nx;
                         dely=(box.yhi-box.ylo)/Ny;
                         delz=(box.zhi-box.zlo)/Nz;
                         centres->vol=delx*dely*delz;
                         x1=box.xlo;y1=box.ylo;z1=box.zlo;
                        for(i=1;i<=Nx;i++){
                               x2=x1+delx;
                            for(j=1;j<=Ny;j++){
                                 y2=y1+dely;       
                              for(k=1;k<=Nz;k++){
                               z2=z1+delz;
                         
                               centres->xc[count]=x1+(x2-x1)*0.5;
                               centres->yc[count]=y1+(y2-y1)*0.5;
                               centres->zc[count]=z1+(z2-z1)*0.5;
                            
                                  z1=z2;
               //printf("%d %lf %lf %lf\n",count,centres->xc[count],centres->yc[count],centres->zc[count]);
                                    count=count+1; 
                          
                                           }
                            
                               z1=box.zlo;
                               y1=y2;
                                         } 
                         
                                x1=x2;
                                y1=box.ylo;
                                      }
                   // printf("%d %lf\n",timestep,delx*dely*delz);
                    return centres;
              
                  }

void free_cubecentres(Centres *centres)
{
free(centres->xc);
free(centres->yc);
free(centres->zc);
free(centres);
}
