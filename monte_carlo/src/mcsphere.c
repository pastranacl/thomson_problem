/****************************************************************************

    mcsphere: Monte-Carlo equilibration of particles in a spherical shell 
            (Thomson's problem) for homogeneous distance between vertexes

    Copyright (C) 2019  Cesar L. Pastrana

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

****************************************************************************/




#include "mcsphere.h"




/* Compile with: gcc mcsphere.c rands.c utils.c -std=c99 -lm --fast-math -fopenmp -O3 -o mcsphere */



int main(int argc, char *argv[])
{
    
    double *rx, *ry, *rz, *E;
    const double R = 1.0;           /* Radious of the Sphere */
    
    if(argc <3) {
        printf("Not sufficient parameters. Usage mcsphere n_particles N_MC_steps\n");
        exit(1);
    }
    
    int N = atoi(argv[1]);                    /* Number of particles */
    int MCSTEPS = atoi(argv[2]);              /* Number of Monte Carlo steps */
    
    rx = (double *)malloc(N*sizeof(double*)); 
    ry = (double *)malloc(N*sizeof(double*));
    rz = (double *)malloc(N*sizeof(double*));
    E =  (double *)malloc(MCSTEPS*sizeof(double*));
    
    printf("Minimization start (N = %d, MC_STEPS = %d: ", N, MCSTEPS);
    sphere(rx, ry, rz, E, R, N,  MCSTEPS);
    
    free(rx); 
    free(ry); 
    free(rz);
    free(E);
    
    return 0;
}


void sphere(double *rx, double *ry, double *rz, double *E, double R, int N, int MC_STEPS)
{
    double *azi, *polar;
    double tx, ty, tz, tpolar, tazi;
    double tE, DE;
    double r;
    double T, dT, kBT, BETA;
    
    int tp;
    long *idum, foornd = -time(NULL);   /* Varibles for randm number generator */   
    idum = &foornd;
    
    char numid[20];    
    char coord_flenme[1000] = COORDS_FILE_NAME;

    FILE *coords_file;
    FILE *energy_file;   
    
    T = INIT_TEMP;
    dT = (INIT_TEMP-1.0)/( (double)MC_STEPS/ANNEAL_TIME );
    kBT=kB*T;
    BETA = 1.0/kBT;
    
    
    /* Initial random positioning of the particles in the surface of a hemisphere */
    azi   = (double *)malloc(N*sizeof(double*)); 
    polar = (double *)malloc(N*sizeof(double*)); 
    for( int i=0; i<N; i++)
    {
        azi[i] = 2*PI*ran1(idum);
        polar[i] = PI*ran1(idum);	
        rx[i]=R*sin(polar[i])*cos(azi[i]);
        ry[i]=R*sin(polar[i])*sin(azi[i]);
        rz[i]=R*cos(polar[i]); 
    }
    
    /* Initial configuration energy */
    tE = 0;
    #pragma omp parallel for reduction(+:tE), private(r)
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            if(i!=j) {
                r = edist( rx[i] - rx[j], ry[i] - ry[j], rz[i] - rz[j] );
                tE += 1.0/r;
            }
        }
    }
    E[0]=tE;
    
    
    /* +++++++++++++++++++++    MAIN MONTE CARLO ROUTINE  +++++++++++++++++++++ */
    for(int n=1; n<MC_STEPS; n++)
    {
        if(n % ANNEAL_TIME==0) { 
            T-=dT; 
            kBT = kB*T; 
            BETA = 1.0/kBT;
            
            // Here you can add if the slope is nearly constant and scape in affirmative case
        }
                
        // Select a random particle and atempt a random move 
        tp = (N - 1)*ran1(idum);
        
        /* Attempted displacement */
        tpolar = polar[tp];
        tazi   = azi[tp];
        tx = rx[tp]; 
        ty = ry[tp]; 
        tz = rz[tp];
        
        polar[tp] += Drad*(ran1(idum)-0.5);
        azi[tp] += Drad*(ran1(idum)-0.5);
        rx[tp] = R*sin(polar[tp])*cos(azi[tp]);
        ry[tp] = R*sin(polar[tp])*sin(azi[tp]);
        rz[tp] = R*cos(polar[tp]);
        

        /* Calculates the energy of the new configuration */
        tE = 0;
        #pragma omp parallel for reduction(+:tE), private(r)
        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
                if(i!=j) {
                    r = edist( rx[i] - rx[j], ry[i] - ry[j], rz[i] - rz[j] );
                    tE += 1.0/r;
                }
            }
        }
        
        /* Metropolis criterion */
        DE = tE - E[n-1];
        if(  ran1(idum) <= exp(-DE*BETA) )
            E[n] = tE;
        else{
            rx[tp] = tx; 
            ry[tp] = ty; 
            rz[tp] = tz;
            azi[tp] = tazi;
            polar[tp] = tpolar;
            E[n] = E[n-1];                   
    }
    
        /* Saves the coordinates of the system after MC optimization */
        /*
        if( n > MCSTEPS - 5000 ){
            strncpy(coord_flenme, COORDS_FILE_NAME, 1000);
            itoa(n, numid);
            strcat(coord_flenme, numid );      
            strcat(coord_flenme, ".dat");
            coords_file = fopen(coord_flenme, "w");  
            for(int i=0; i<N; i++)
                fprintf(coords_file, "%.10g\t%.10g\t%.10g\n", rx[i], ry[i], rz[i]);
            fclose(coords_file); 
        }*/
    }
    /* ---------------------------------------------------------------------------- */

    free(polar);
    free(azi);
    
    
    // DATA SAVING 
    /* Writes the energy of the system as a function of the MC steps */
    energy_file = fopen(MC_ENRGY_FILE_NAME, "w");
    for(int n=0; n<MC_STEPS; n+=200)
        fprintf(energy_file, "%d\t%.10g\n", n, E[n]); 
    fclose(energy_file);

    /* Writes the eq. configuration  */    
    coords_file = fopen(COORDS_FILE_NAME, "w");  
    for(int i=0; i<N; i++)
        fprintf(coords_file, "%.10g\t%.10g\t%.10g\n", rx[i], ry[i], rz[i]);

    printf("Done! Final Temp. T = %f K\n", T);
    fclose(coords_file);

    
}


