#include "stdheader.h"
int getnatoms(FILE *fp) {
int i,natoms;
size_t nbytes = 100; 
char *my_string,str1[100], str2[100];
my_string = (char *) malloc(nbytes + 1);
		for (i = 0; i < 5; i++) {
			my_string = (char *) malloc(nbytes + 1);
			getline(&my_string, &nbytes, fp);
			sscanf(my_string, "%s %s", str1, str2);
			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "NUMBER") == 0) {
				getline(&my_string, &nbytes, fp);
				sscanf(my_string, "%d", &natoms);
			}
		}
                free(my_string);
		return natoms;
	}


