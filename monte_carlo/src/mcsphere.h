#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "rands.h"
#include "utils.h"

#define PI 	3.141592653589793238462643383279502884197169
#define kB 1.38064852*1e-23

#define Drad 1.0*PI/180.0  
#define INIT_TEMP 350
#define ANNEAL_TIME 50

#define MC_ENRGY_FILE_NAME  "energy.dat"
#define COORDS_FILE_NAME    "rsphere_eq.dat"


void sphere(double *rx, double *ry, double *rz, double *E, double R, int N, int MC_STEPS);


inline double edist(double x, double y, double z)
{
    return x*x + y*y + z*z; 
}
