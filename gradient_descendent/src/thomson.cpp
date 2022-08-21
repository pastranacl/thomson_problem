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
#include "thomson.hpp"



/**************************************************************************
*
* Particles: Particles constructor
* 
* Input:        N = Number of particles
*               R = Radius of the spherical surface
**************************************************************************/
Particles::Particles(int N, float R)
{ 
    rSpherical.resize(N*2); 
    rCartesian.resize(N*2); 
    this->R = R;
    this->N = N; // Equivalent to Particles::N = N
    lb = bondLength();
    setInitialConfiguration();
}


/**************************************************************************
*
* setInitialConfiguration: Minimise the repulsion energy between points by 
* using the Barzilai-Borwain's gradient descendent implementation.
* 
* Input:        maxits = Maximum number of iterations
*               gtol   = Tolerance (norm of the gradient vector) to stop 
*                        the evaluation
* 
*--------------------------------------------------------------------------
*
*  Ouput: The gradient along the phi and theta directions by considering
*         a Coulombic repulstion between particles
*
**************************************************************************/
void Particles::setInitialConfiguration()
{
    int np=0;
    double l0Char = FACT_L0_INIT*lb;
    float polar0 = 0.5*lb;
    srand (time(NULL));
    rSpherical(0) = randInterval(0.0, 2.0*PI);
    rSpherical(1) = randInterval(polar0, PI-polar0);
    
    np=1;
    do {
        
        rSpherical(np*2)   = randInterval(0.0, 2.0*PI);
        rSpherical(np*2+1) = acos(randInterval(-1.0, 2.0)); 
        
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

/**************************************************************************
*
* bondLength: Characteristic length between particles in the unit sphere
* 
*--------------------------------------------------------------------------
*
*  Ouput: Characteristic length between particles in the unit sphere
*
**************************************************************************/
double Particles::bondLength()
{
    return sqrt(8.0*PI/((N-2)*sqrt(3.0)) );
}


/**************************************************************************
*
* minimise: Minimise the repulsion energy between points by using the 
* Barzilai-Borwain's gradient descendent implementation.
* 
* Input:        maxits = Maximum number of iterations
*               gtol   = Tolerance (norm of the gradient vector) to stop 
*                        the evaluation
* 
*--------------------------------------------------------------------------
*
*  Ouput: The gradient along the phi and theta directions by considering
*         a Coulombic repulstion between particles
*
**************************************************************************/
void Particles::minimise(int maxits, double gtol, bool saveIts)
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
    
    Export dump(this, "./output/");
    
    
    // Iterative minimisations
    for(int i=0; i<maxits; i++) {
        
        // Calculate gradient 
        calcGradient(grad);
        if(vecNanOrInf(grad)) break;
            
        // Calculate damping factor
        Dr = rSpherical - rSpherical_pre;
        Dgrad = grad - grad_pre;
        
        float gamma = abs(Dr.dot(Dgrad))/Dgrad.squaredNorm();
       
        rSpherical_pre = rSpherical;
        grad_pre = grad;
        
        // Move
        rSpherical -= gamma*grad;
        
        if(vecNanOrInf(rSpherical)) {
            rSpherical = rSpherical_pre;
            break;
        }

        // Scape
        if(grad.norm()<gtol)
            break;

        
        // Save snapshots if applicable
        if(saveIts) 
            dump.exportXYZ(dump.getItfileName(i, "sphere"));
    }
    
}


/**************************************************************************
*
* calcGradient: Check if the input vector has any infinite value
*--------------------------------------------------------------------------
*
*  Ouput: The gradient along the phi and theta directions by considering
*         a Coulombic repulstion between particles
*
**************************************************************************/
void Particles::calcGradient(VectorXf& grad)
{
    grad.setZero();
    #pragma omp parallel for 
    for(int i=0; i<N; i++) {
        
        float phi1, theta1;
        phi1 = rSpherical(i*2);
        theta1 = rSpherical(i*2+1);
        for(int j=0; j<N; j++) {
            if(i==j) continue;
            
            float phi2, theta2;
            phi2 = rSpherical(j*2);
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


/**************************************************************************
*
* gradInf: Check if the input vector has any infinite value
*--------------------------------------------------------------------------
*
*  Ouput: false if no Inf values found
*
**************************************************************************/
bool Particles::vecNanOrInf(VectorXf& v)
{
    for(int i=0;i<N; i++){
        if(isinf(v(i*2)) || isnan(v(i*2)))
            return true;
        if(isinf(v(i*2+1)) || isnan(v(i*2+1)))
            return true;
    }
    
    return false;
}

/**************************************************************************
*
* getCartesian: Converts the spherical coordiantes to cartesian (privat)
*--------------------------------------------------------------------------
*
*  Ouput: Updates the protected function rCartesian and returns the value
*
**************************************************************************/
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



/**************************************************************************
*
* Export constructor
* 
* Input: particles = Class to save the Data
*        outPath   = Path to save the data
* 
**************************************************************************/
Export::Export(Particles *particles, std::string outPath)
{ 
    this->particles = particles;
    this->outPath = outPath;
}


/**************************************************************************
*
* exportVTU: Exports the data to an VTK XML UnstructuredGrid, such that 
*            the point particles can be visualised with ParaView
*
*  Input: Name of file
*
**************************************************************************/
void Export::exportVTU(std::string fileName)
{
    std::ofstream fid;
    fid.open(outPath + fileName + ".vtu");
    
    double d = 2.0*particles->R*particles->lb;

    fid << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    fid << "  <UnstructuredGrid>" << std::endl;
    fid << "    <Piece NumberOfPoints=\"" << particles->N << "\" NumberOfCells=\"0\">" << std::endl;
    fid << "      <Points>" << std::endl;
    fid << "        <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    fid << "        ";
    for(int i=0; i<particles->N; i++) {
        float phi   = particles->rSpherical(2*i), theta = particles->rSpherical(2*i+1);
        fid << particles->R*sin(theta)*cos(phi) << " "
            << particles->R*sin(theta)*sin(phi) << " " 
            << particles->R*cos(theta) << " ";
    }
    fid << std::endl;
    fid << "        </DataArray>" << std::endl;
    fid << "      </Points>" << std::endl;
    fid << "      <PointData>" << std::endl;
    fid << "        <DataArray type=\"Float32\" Name=\"Diameter\" format=\"ascii\">" << std::endl; 
    fid << "        ";
    for(int i=0; i<particles->N; i++)     fid << d << " "; fid << std::endl;
    fid << "        </DataArray>" << std::endl;
    fid << "      </PointData>" << std::endl;
    fid << "      <Cells>" << std::endl;
    fid <<"         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    fid <<"         </DataArray>" << std::endl;
    fid <<"         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    fid <<"         </DataArray>" << std::endl;
    fid <<"         <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    fid <<"         </DataArray>" << std::endl;
    fid <<"       </Cells>" << std::endl;
    fid << "    </Piece>" << std::endl;
    fid << "  </UnstructuredGrid>" << std::endl;
    fid << "</VTKFile>";
    fid.close();
}

/**************************************************************************
*
* exportXYZ: Exports the data to XYZ file. Simple file, allow simple 
*            visualisation with Ovito
*
*  Input:     fileName =  Name of the file
*
**************************************************************************/
void Export::exportXYZ(std::string fileName)
{
    std::ofstream fid;
    fid.open(outPath + fileName + ".xyz");
    fid << particles->N << std::endl << std::endl;
    for(int i=0; i<particles->N; i++) {
        float phi   = particles->rSpherical(2*i),  theta = particles->rSpherical(2*i+1);
        fid << particles->R*sin(theta)*cos(phi) << "\t"
            << particles->R*sin(theta)*sin(phi) << "\t" 
            << particles->R*cos(theta) << std::endl;
    }
    fid.close();
}


/**************************************************************************
*
* exportDAT: Exports the data to XYZ file. Simplest file
*
*  Input:   fileName =  Name of the file
*
**************************************************************************/
void Export::exportDAT(std::string fileName)
{
    std::ofstream fid;
    fid.open(outPath + fileName + ".dat");
    for(int i=0; i<particles->N; i++) {
        float phi   = particles->rSpherical(2*i),  theta = particles->rSpherical(2*i+1);
        fid << particles->R*sin(theta)*cos(phi) << "\t"
            << particles->R*sin(theta)*sin(phi) << "\t" 
            << particles->R*cos(theta) << std::endl;
    }
    fid.close();
}




/**************************************************************************
*
* getItfileName: Genrates full filename for the minimisation iterations
*
*  Input:        it = Integer preceeding the main name
*                fileName = Main name
* 
*--------------------------------------------------------------------------
* Output:        fullName = String including the full path and the it, e.g.
*                           ./path/1_file_name
* 
**************************************************************************/
std::string Export::getItfileName(int it, std::string fileName)
{
    std::stringstream ss;
    ss << it << "_" << fileName;
    std::string fullName = ss.str();
    return fullName;    
}
