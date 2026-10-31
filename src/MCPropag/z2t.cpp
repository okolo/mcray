#include "Cosmology.h"
#include "Units.h"
#include <iostream>
using namespace MCPropag;

int main(int argc, char** argv) {
    //Test::PrecisionTests::SetPrecision(128);
    int result = 0;
    try
    {
        double z=-1.;
        if(argc<2 || sscanf(argv[1],"%lg",&z)!=1 || z<0) {
            std::cerr << "Calculate light travel time (in Mpc) from remote object\n"
            << "usage: " << argv[0] << "<z>" << std::endl;
            return 1;
        }
        if(!cosmology.IsInitialized())
            cosmology.Init(z+10.);
        double t = (cosmology.getAgeOfUniverse()-cosmology.z2t(z));
        std::cout << t/units.Mpc << std::endl;
    }
    catch(Exception* ex)
    {
        std::cerr << ex->Message() << std::endl;
        result = 1;
    }
    return result;
}