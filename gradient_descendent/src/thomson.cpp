#include "thomson.hpp"



Particles::Particles(int N, float R)
{ 
    rSpherical.resize(N*2); 
    rCartesian.resize(N*2); 
    this->R = R;
    this->N = N; // Equivalent to Particles::N = N
    setInitialConfiguration();
}


void Particles::setInitialConfiguration()
{
    int np=0;
    double l0Char = FACT_L0_INIT*bondLength();
    
    srand (time(NULL));
    rSpherical(0) = rand_zero2Val(2.0*PI);
    rSpherical(1) = rand_zero2Val(PI);
    
    np=1;
    do {
        
        rSpherical(np*2)   = rand_zero2Val(2.0*PI);
        rSpherical(np*2+1) = rand_zero2Val(PI);
        
        float td;
        for(int i=0; i<np; i++) {
            td = acos(dotSpherical(rSpherical(i*2),  rSpherical(i*2+1), 
                                   rSpherical(np*2), rSpherical(np*2+1)));
            if(td<l0Char)
                break;
        }
        if(td>l0Char)
            np++;
        
        
    } while(np<N);
}

double Particles::bondLength()
{
    return sqrt(8.0*PI/((N-2)*sqrt(3.0)) );
}


void Particles::minimise(int maxits, double gtol)
{
    
    VectorXf rSpherical_pre(N*2);
    VectorXf grad(N*2);
    VectorXf grad_pre(N*2);
    VectorXf Dr(N*2);
    VectorXf Dgrad(N*2);
    
  
    // Initial displacement
    calcGradient(grad);
    rSpherical_pre = rSpherical;
    grad_pre = grad;
    rSpherical -= GAMMA_0*grad;
    
        
    // Iterative minimisations
    for(int i=0; i<maxits; i++) {
        float gamma;
        // Calculate gradient 
        calcGradient(grad);
        if(gradInf(grad)) break;
            
        // Calculate damping factor
        Dr = rSpherical - rSpherical_pre;
        Dgrad = grad - grad_pre;
        gamma = abs(Dr.dot(Dgrad))/Dgrad.squaredNorm();
       
        rSpherical_pre = rSpherical;
        grad_pre = grad;
        
        // Move
        rSpherical -= gamma*grad;
        

        // Scape
        if(grad.norm()<gtol)
            break;

        
        //
        std::stringstream ss;
        ss << "./output/" << i << "_file.xyz";
        std::string fname = ss.str();
        exportXYZfile(fname);
    }
    
}


void Particles::calcGradient(VectorXf& grad)
{
    grad.setZero();
    #pragma omp parallel for 
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            if(i==j) continue;
            
            float phi1, phi2, theta1, theta2;
            phi1 = rSpherical(i*2);
            phi2 = rSpherical(j*2);
            theta1 = rSpherical(i*2+1);
            theta2 = rSpherical(j*2+1);
            
            float dotp, alpha, pref;
            float sinTh2cosPh2, sinTh2sinPh2;
            dotp = dotSpherical(phi1, theta1, phi2, theta2);
            alpha = acos(dotp);
            pref = 1.0/(alpha*alpha*sqrt(1-dotp*dotp));
            
            sinTh2cosPh2 = sin(theta2)*cos(phi2);
            sinTh2sinPh2 = sin(theta2)*sin(phi2);
            
            #pragma omp critical
            grad(i*2+1)  +=  pref*(cos(theta1)*cos(phi1)*sinTh2cosPh2 
                                 + cos(theta1)*sin(phi1)*sinTh2sinPh2 
                                 - sin(theta1)*cos(theta2));
            
            #pragma omp critical
            grad(i*2) += pref*( - sin(theta1)*sin(phi1)*sinTh2cosPh2 
                                + sin(theta1)*cos(phi1)*sinTh2sinPh2);
        }
    }
}

bool Particles::gradInf(VectorXf& grad)
{
    for(int i=0;i<N; i++){
        if(isinf(grad(i*2)) )
            return true;
        if(isinf(grad(i*2+1)) )
            return true;
    }
    
    return false;
}


VectorXf Particles::getCartesian()
{
    for(int i=0;i<N; i++)
    {
        double phi   = rSpherical(2*i),
               theta = rSpherical(2*i+1);
        rCartesian(3*i)   = sin(theta)*cos(phi);
        rCartesian(3*i+1) = sin(theta)*sin(phi);
        rCartesian(3*i+2) = cos(theta);
    }
    
    return R*rCartesian;
}





void  Export::exportXYZfile(std::string fileName)
{
    std::ofstream fid;
    fid.open(fileName);
    fid << N << std::endl << std::endl;
    for(int i=0; i<N; i++) {
        double phi   = rSpherical(2*i),
               theta = rSpherical(2*i+1);
        fid << R*sin(theta)*cos(phi) << "\t"
                    << R*sin(theta)*sin(phi) << "\t" 
                    << R*cos(theta) << std::endl;
    }
    fid.close();
}
