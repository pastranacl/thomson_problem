#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define ABS(a)  (a) < 0 ? -(a) : (a)
#define DM  3

#define     PI           3.14159265358979323846

#define     INIT_STEP_BACK_TRACK    1.0e-5          // Inital step for backtracking
#define     MAX_ITS_LINE_SEARCH     100             //

#define     MAX_STEPS_GRAD_DESC    5

#define     K_SPHERE     0.25                         // Strength of the spherical potential
#define     KC           0.0                       // Strength of the Coulombic potential
#define     R_CUT        1.0                        // Cutoff distance for the calculation of the particle-particle interaction

#define     FNAME_COORDS_MINIM  "minim_coords.dat"



void sphere_init_coords(double *x, int N);
void backtrack(double *x, double *g, int N);
double energy(double *x, int N);

void grad(double *r, double *f, int N);
void thomson_gdbb(double *x, int N);

void daddvec(double *v, double *u, double a, double *w, int N);
double ddot(double *v, double *u, int N);
double dsetzerovec(double *v, int N);
