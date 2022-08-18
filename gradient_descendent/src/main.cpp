#include <iostream>
#include "thomson.hpp"



int main()
{
    
    int N;
    double R;
    std::cout << "Set number of particles: "; std::cin >> N;
    std::cout << "Set sphere radius: "; std::cin >> R;
    
    Particles surface(N, R);
    surface.minimise(100, 1.0e-6);
    //surface.exportXYZfile("test.xyz");
    
    return 0;
}

 