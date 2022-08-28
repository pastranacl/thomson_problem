/***************************************************************************
    thomson: Generates a homogeneous distribution of particles on a spherical 
             surface using gradient minimisation

    Copyright (C) Cesar L. Pastrana, 2022

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

***************************************************************************/
#include <iostream>
#include <fstream>     /* ofstream */
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <Eigen/Dense>
#include <Eigen/Core>   /*for seq*/

#ifndef THOMSON
#define THOMSON

#define PI              3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define FACT_L0_INIT    1.0e-1
#define GAMMA_0         1.0e-6


typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorXf;

 
class Particles {
    public:
        Particles(int N, float R, std::string initProcedure = "random");
        VectorXf rSpherical; // if we knew the size in advance: rSpherical(n)
        VectorXf rCartesian;
        float R;
        int N;
        float lb;
       
        void minimise(int maxits, double tol, bool saveIts = true);
        void calcGradient(VectorXf& grad);
        VectorXf getCartesian();
    
    protected:
        double bondLength();
        
    private:
        void setInitialConfiguration_Random();
        void setInitialConfiguration_Fibonacci();
        bool vecNanOrInf(VectorXf& v);
        inline float dotSpherical(float phi1, float theta1, 
                                  float phi2, float theta2);
        inline float randInterval(float min_val, float max_val);
};




class Export {

    public:
        Export(Particles *particles, std::string outPath);
        Particles *particles;
        std::string outPath;
        void exportXYZ(std::string fileName);
        void exportVTU(std::string fileName);
        void exportDAT(std::string fileName);
        std::string getItfileName(int it, std::string fileName);
        
    private:
        std::string fileNameIt(int id);
};




class Import {
    public:
        Import(Particles *particles, std::string inPath);        
        Particles *particles;
        std::string inPath;
        void importXYZ(std::string fileName);
        void importDAT(std::string fileName);
};
 




/**************************************************************************
*
* dotSpherical: Calculate the dot product between a couple points r and u
*               in a unit spherical surface
* 
* Input:        phi_r   = Azimuth angle particle r
*               theta_r = Polar angle first particle r
*               phi_u   = Azimuth angle particle u
*               theta_u = Polar angle first particle u
* 
*--------------------------------------------------------------------------
*
*  Ouput: Dot product
*
**************************************************************************/
inline float Particles::dotSpherical(float phi1, float theta1, 
                                     float phi2, float theta2)
{
    float prodSinTheta;
    prodSinTheta = sin(theta1)*sin(theta2);
    return prodSinTheta*cos(phi1)*cos(phi2) + prodSinTheta*sin(phi1)*sin(phi2) + cos(theta1)*cos(theta2);
}


/**************************************************************************
*
* rand_zero2Val: Generates random numbers between min_val and max_val
* 
* Input:        max_val   =  Maximum value of the random numbers interval
* 
*--------------------------------------------------------------------------
*
*  Ouput: Random number in the specified interval
*
**************************************************************************/
inline float Particles::randInterval(float min_val, float max_val)
{
    return min_val + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(max_val)));
} 


#endif
