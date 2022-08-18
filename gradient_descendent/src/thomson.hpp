#include <fstream>     /* ofstream */
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <Eigen/Dense>

#define PI              3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define FACT_L0_INIT    5.0e-1
#define GAMMA_0         1.0e-6



typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorXf;

 
class Particles {
    public:
        VectorXf rSpherical; // if we knew the size in advance: rSpherical(n)
        float R;
        int N;
        
        Particles(int N, float R);
        void minimise(int maxits, double tol);
        void calcGradient(VectorXf& grad);
        void exportXYZfile(std::string fileName);
        VectorXf getCartesian();
        
    private:
        VectorXf rCartesian;
        void setInitialConfiguration();
        double bondLength();
        bool gradInf(VectorXf& grad);
        inline float dotSpherical(float phi1, float theta1, 
                                  float phi2, float theta2);
        inline float rand_zero2Val(float max_val);
        
};


class Export: public Particles {

    public:
    private:
        string fileNameIt(int id);
    // Export XYZ 
    // COMBINE STRING
    // EXPORT OTHER FORMAT
};


inline float Particles::dotSpherical(float phi1, float theta1, 
                                     float phi2, float theta2)
{
    float prodSinTheta;
    prodSinTheta = sin(theta1)*sin(theta2);
    return prodSinTheta*cos(phi1)*cos(phi2) + prodSinTheta*sin(phi1)*sin(phi2) + cos(theta1)*cos(theta2);
}


inline float Particles::rand_zero2Val(float max_val)
{
    return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(max_val)));
} 
