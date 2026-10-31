// Angle inclusion selects modern Jet.h, not the extracted reference header.
#include <Jet.h>
#include "fork_api.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("Usage: fork-vmf-weight-probe POINTS");
        std::ifstream input(argv[1]);
        if(!input) throw std::runtime_error("Missing weight points");
        std::cout << std::setprecision(21) << std::boolalpha;
        std::string group;
        double theta,phi,kappa;
        long double direction[3];
        while(input >> group >> theta >> phi >> kappa >> direction[0] >> direction[1] >> direction[2]) {
            const auto old=ForkVmfWeight(theta,phi,direction,kappa);
            mcray::coord_type modernDirection[3]{direction[0],direction[1],direction[2]};
            const auto modern=mcray::vMF_weight(theta,phi,modernDirection,kappa);
            std::cout << "{\"group\":\"" << group << "\",\"theta\":" << theta << ",\"phi\":" << phi
                      << ",\"kappa\":" << kappa << ",\"fork_finite\":" << std::isfinite(old)
                      << ",\"fork\":";
            if(std::isfinite(old)) std::cout << old; else std::cout << "null";
            std::cout << ",\"modern\":" << modern << ",\"relative_difference\":";
            if(std::isfinite(old)) std::cout << std::abs(old-modern)/std::max({std::abs(old),std::abs(modern),1e-300});
            else std::cout << "null";
            std::cout << "}\n";
        }
        if(!input.eof()) throw std::runtime_error("Malformed weight points");
        return 0;
    } catch(Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error); std::cerr << error->Message() << '\n';
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
