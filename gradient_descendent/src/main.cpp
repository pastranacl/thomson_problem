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
    Particles surface(N, R);
    surface.minimise(100, 1.0e-6);
    Export dump(&surface, "./");
    dump.exportDAT("rsphere_eq");
    
    std::cout << "Done!" << std::endl;
    return 0;
}

 
