#include <iostream>

#include "thomson.hpp"





int main(int argc, char *argv[])
{
     
    if(argc <3) {
        printf("Not sufficient parameters. Usage thomson n_particles radius\n");
        exit(1);
    }
 
    int N = atoi(argv[1]);                  /* Number of particles */
    double R = atof(argv[2]);              /* Radius */
    
    std::cout << "Minimising... ";
    
    Particles surface(N, R, "random");
    
    //Import importer(&surface, "./");
    //importer.importXYZ("init_coords.xyz");
    
    surface.minimise(2000, 1.0e-8);
    
    Export dump(&surface, "./");
    dump.exportDAT("rsphere_eq");
    
    std::cout << "Done!" << std::endl;
    
    return 0;
}

 
