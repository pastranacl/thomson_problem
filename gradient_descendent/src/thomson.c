#include "thomson.h"
#include "rands.h"




int main()
{
    
    
    //Allocate x and position inital coordinates randomly on spherical shell
    int N=5;
    int n_eff = DM*N;
    double *x;
    x = (double *)malloc(n_eff*sizeof(double));
    sphere_init_coords(x, N);
    
    
    // Minimisation
    double time_exec;
    clock_t tic, toc; 
    
    tic = clock(); 
    printf("Minimisation... ");
    
    thomson_gdbb(x, N);
    toc = clock();
    time_exec = (double)(toc - tic)/CLOCKS_PER_SEC;       
    printf("Done (elapsed time = %f seconds)\n", time_exec);
    
    // Save data
    FILE *fidcords=fopen(FNAME_COORDS_MINIM, "w+");
    for(int i=0;i<N; i++)
        fprintf(fidcords, "%f\t%f\t%f\n", x[i*DM], x[i*DM+1], x[i*DM+2]);
    
    return 0;
}




void thomson_gdbb(double *x, int N)
{

    int n_eff = DM*N;
    double *df, *dx;
    double *x_pre;
    double *f, *f_pre;
    
    df = (double *)malloc(n_eff*sizeof(double));
    dx = (double *)malloc(n_eff*sizeof(double));
    x_pre = (double *)malloc(n_eff*sizeof(double));
    f     = (double *)malloc(n_eff*sizeof(double));
    f_pre = (double *)malloc(n_eff*sizeof(double));
    
    
    // First step is a simple backtrack line-search
    memcpy(x_pre, x, n_eff*sizeof(double));
    grad(x_pre, f_pre, N);
    backtrack(x, f_pre, N);
    

    // Barzilai-Borwain Gradient-descendent 
    for(int t=0; t<MAX_STEPS_GRAD_DESC*N; t++) {

        grad(x, f, N);
        
       
        double proj_dxdf, dfnorm;
        double gamma; 
        daddvec(f, f_pre, -1.0, df, n_eff);
        daddvec(x, x_pre, -1.0, dx, n_eff);
        proj_dxdf = ddot(dx, df, n_eff);             
        dfnorm = ddot(df, df, n_eff);
        
        if(dfnorm==0.0) break;
        gamma = ABS(proj_dxdf/dfnorm);
        
            
        memcpy(x_pre, x, n_eff*sizeof(double *));    
        memcpy(f_pre, f, n_eff*sizeof(double *));    
        

       

        daddvec(x_pre, f, -gamma, x, n_eff);   
         for(int i=0; i<N; i++)
            printf("%f\t%f\t%f\n", x[i*DM], x[i*DM+1], x[i*DM + 2]);
         printf("\n\n");
       //break;
    }
    
   
    free(f);
    free(f_pre);
    free(x_pre);
}




/*********************************************************/
/*                          grad                         */
/* Calculate the gradient as a result of reppulsion and  */
/* confinement to the unit sphere                        */
/*********************************************************/
void grad(double *r, double *f, int N)
{
    
    dsetzerovec(f, N);
    
    for(int i=0; i<N; i++) {
        
        double xi, yi, zi;
        xi = r[i*DM];
        yi = r[i*DM+1];
        zi = r[i*DM+2];
        
        for(int j=i+1; j<N; j++){
            
            double xj, yj, zj;
            double dx, dy,dz, d, id, id3;
            double fc;
            xj = r[j*DM];
            yj = r[j*DM+1];
            zj = r[j*DM+2];
            dx = xi - xj,
            dy = yi - yj,
            dz = zi - zj,
            d = dx*dx + dy*dy + dz*dz;
            if (i==j || d>0.5) continue;
            
            id = 1/d;
            id3 = id*id*id;
            fc = KC*id3;
            
            f[i*DM]   -= fc*dx;
            f[i*DM+1] -= fc*dy;
            f[i*DM+2] -= fc*dz;
            
            
            f[j*DM]   += fc*dx;
            f[j*DM+1] += fc*dy;
            f[j*DM+2] += fc*dz;
        }
        
        
        // Radius constrain
        double ri = sqrt(xi*xi + yi*yi + zi*zi);
        double f_sp = K_SPHERE*(ri - 1.0);    
        
        f[i*DM]   += f_sp*xi/ri;
        f[i*DM+1] += f_sp*yi/ri;
        f[i*DM+2] += f_sp*zi/ri;
        
    }
    
     dsetzerovec(f, N);
    
}



double energy(double *x, int N)
{    
    double E=0;
    
    for(int i=0; i<N; i++) {
        double xi, yi, zi;
        xi = x[i*DM];
        yi = x[i*DM+1];
        zi = x[i*DM+2];
        for(int j=i+1; j<N; j++) {
            double xj, yj, zj;
            double dx, dy,dz, d, id, id3;
            double fc;
            xj = x[j*DM];
            yj = x[j*DM+1];
            zj = x[j*DM+2];
            dx = xi - xj,
            dy = yi - yj,
            dz = zi - zj,
            d = dx*dx + dy*dy + dz*dz;
            E += 1/(2*d);
        }
        
        // Energy confinement
        double r =sqrt(xi*xi + yi*yi + zi*zi);
        double dr = (r-1.0);
        E += K_SPHERE*dr*dr/2;
    }
    
    return E;
}


/***************************************************************/
/*                          backtrack                          */
/* Back-track line search: find the factor alpha minimizing    */
/* the energy.                                                 */
/***************************************************************/
void backtrack(double *x, double *g, int N)
{
    
    double E_pre, E;
    double alpha = INIT_STEP_BACK_TRACK;
    int its;
    
    double *tx;
    tx = (double *)malloc(N*DM*sizeof(double));
    memcpy(tx, x, N*DM*sizeof(double));
    
    
    E = energy(tx, N);
    its = 0;
    do
    {
        E_pre = E;
        daddvec(x, g, alpha, tx, DM*N);
        
        E = energy(tx, N);
        alpha *= 1.01;
        its++;
                
    } while(E_pre > E && its<MAX_ITS_LINE_SEARCH);
    
    // We cannot use daddvec over x and x (racing condition)
    #pragma omp parallel for
    for(int i=0; i<N*DM; i++)
        x[i] += alpha*g[i];
}

/***************************************************************/
/*                      sphere_init_coords                     */
/* Inital coordinates randomly initialized                     */
/***************************************************************/
void sphere_init_coords(double *x, int N)
{
    
    long *idum, foornd;                 /* Varibles for randm number generator */  
    foornd = -time(NULL);
    idum = &foornd;
    
    for(int i=0; i<N; i++) {
        double phi, theta;
        phi = 2*PI*ran1(idum);
        theta = PI*ran1(idum);
        x[i*DM]   = cos(theta)*sin(phi);
        x[i*DM+1] = sin(theta)*sin(phi);
        x[i*DM+2] = cos(phi);
    }
}




/***************************************************************/
/*                          daddvec                            */
/* Sums two vectors v,u in linear combination w =  v + a*u     */
/***************************************************************/
void daddvec(double *v, double *u, double a, double *w, int N)
{
    #pragma omp parallel for
    for(int i=0; i<N; i++) 
        w[i] = v[i] + a*u[i];
}



/***************************************************************/
/*                          ddot                               */
/* Dot product between vecors u abd v, (u·v)                   */
/***************************************************************/
double ddot(double *v, double *u, int N)
{
    int dp=0;
    
    #pragma omp parallel for reduction(+:dp)
    for(int i=0; i<N; i++)
        dp += v[i]*u[i];
    return dp;
}

/***************************************************************/
/*                         dsetzerovec                         */
/* Set all the entries in the vector to 0                      */
/***************************************************************/
double dsetzerovec(double *v, int N)
{
     #pragma omp parallel for
    for(int i=0; i<N; i++)
        v[i] = 0;
}
 
