#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

double erand48 (unsigned short xsubi[3]);

void main(int argc, char *argv[])

{	int nsites,d,i;
	unsigned short X[3];
	if (!argv[1] || !argv[2]) exit(3);
	sscanf(argv[1], "%d", &nsites);
	sscanf(argv[2], "%d", &d);
	X[1] = time(0);
	while(nsites>0){
		for (i=0;i<d;i++) printf("%6.0f ", floor(1e6*erand48(X)));
/*		for (i=0;i<d;i++) printf("%G ", erand48(X));	*/
		printf("\n");
		nsites--;
	}
}
