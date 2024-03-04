#include "stdheader.h"
int gettimestep(FILE *fp)
{
int timestep;
	size_t nbytes = 100;
	char *my_string=NULL;
getline(&my_string, &nbytes, fp);
sscanf(my_string, "%d", &timestep);
//getline(my_string,nbytes, fp);
//sscanf(my_string, "%d", &timestep);

return timestep;
}

