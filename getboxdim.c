#include "stdheader.h"
#include "getboxdim.h"
Box boxdim(FILE *fp)
{
                         
size_t nbytes = 100;
	char *my_string=NULL;
                       Box box;
                       getline(&my_string, &nbytes, fp);
			sscanf(my_string, "%lf %lf", &box.xlo, &box.xhi);
			//printf ("%lf %lf\n",box->xlo,box->xhi);
			getline(&my_string, &nbytes, fp);
			sscanf(my_string, "%lf %lf", &box.ylo, &box.yhi);
			//printf ("%lf %lf\n",box->ylo,box->yhi);
			getline(&my_string, &nbytes, fp);
			sscanf(my_string, "%lf %lf", &box.zlo, &box.zhi);
			//printf ("%lf %lf\n",box->zlo,box->zhi);
                        
                        return box;

} 



