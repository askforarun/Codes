#include "stdheader.h"
#include "Atom.h"
#include "getboxdim.h"
Atom *atominfo(FILE *fp1,Box box,Ljcoeff *lj, int natoms) {
			
                                size_t nbytes = 100;
	                        char *my_string=NULL;
                                int i;
                         	Atom *atom=(Atom *) malloc(sizeof(Atom));
                                atom->id = (int *) malloc(sizeof(int) *natoms);
                                atom->type =(int *) malloc(sizeof(int) * natoms);
				atom->x = (double *) malloc(sizeof(double) * natoms);
				atom->y = (double *) malloc(sizeof(double) * natoms);
				atom->z = (double *) malloc(sizeof(double) * natoms);
				atom->q = (double *) malloc(sizeof(double) * natoms);
                                atom->eps=(double *) malloc(sizeof(double) * natoms);
                                atom->sigma=(double *) malloc(sizeof(double) * natoms);
                                atom->natoms=natoms;
                                	
				for (i = 0; i < natoms; i++) {
                                        atom->eps[i]=lj->eps[i];
                                        atom->sigma[i]=lj->sigma[i]; 
      					getline(&my_string, &nbytes, fp1);
                                        
					sscanf(my_string, "%d %d %lf %lf %lf %lf\n", &atom->id[i],&atom->type[i],&atom->x[i],
							&atom->y[i], &atom->z[i], &atom->q[i]);

                                       if((atom->x[i] - box.xlo) <  0) {atom->x[i]= atom->x[i]+(box.xhi-box.xlo);}//add box length
                                       if((atom->x[i] - box.xhi) >= 0) {atom->x[i]= atom->x[i]-(box.xhi-box.xlo);}//subtract box length
                                       if((atom->y[i] - box.ylo) <  0) {atom->y[i]= atom->y[i]+(box.yhi-box.ylo);}
                                       if((atom->y[i] - box.yhi) >= 0) {atom->y[i]= atom->y[i]-(box.yhi-box.ylo);} 
                                       if((atom->z[i] - box.zlo) <  0) {atom->z[i]= atom->z[i]+(box.zhi-box.zlo);}
                                       if((atom->z[i] - box.zhi) >= 0) {atom->z[i]= atom->z[i]-(box.zhi-box.zlo);}
                                      
    
                                      
                                    }
                                 
                                       
                                      
                                            
                                     
                                         
                                       
                                      
                                         
                                        
                                        

                                          
                                     
                                        
                             
                                      
                                        
                                     
                              
                                         
                                          
                                                 
//                                              fprintf(fp3,"%d %lf %lf %lf\n",atom->id[i],atom->x[i],atom->y[i],atom->z[i]);
                                          // fprintf("%d %lf %lf %lf %lf\n",atom->id[i],atom->x[i],atom->y[i],atom->z[i],(0.5)*atom->sigma[i]);
//                       				printf("%d %lf %lf %lf %lf\n",atom->type[i],atom->x[i],atom->y[i],atom->z[i],atom->sigma[i]);
                                               
                                       
                                              
                                             
                                            
			                           
         
				return atom;
                              
			
}


		
void free_atom(Atom *atom)
{
free(atom->x);
free(atom->type);
free(atom->id);
free(atom->y);
free(atom->z);
free(atom->q);
free(atom->eps);
free(atom->sigma);
free(atom);
}
