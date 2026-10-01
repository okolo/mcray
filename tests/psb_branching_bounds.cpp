#include "Nucleus.h"
#include <climits>
#include <iostream>
#include <stdexcept>

static void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    using namespace Interactions;
    try {
        const double coefficients[] = {0.1, 0.3, 0.6};
        for (int count = 1; count <= 3; ++count) {
            const TBranching branching = {count, coefficients, -1., -1.};
            const TIsotope isotope = {56, 26,
                {{10., 20., 1., 5., 5.}, {10., 20., 1., 5., 5.}},
                1., &branching};
            CPhotoDisintegrationCanal multiple(isotope, EMultiple);
            const int invalid[] = {INT_MIN, -1, 0, count + 1, INT_MAX};
            for (int delta : invalid)
                check(multiple.getRate(delta) == 0., "invalid deltaA must return zero");
            for (int delta = 1; delta <= count; ++delta)
                check(multiple.getRate(delta) == coefficients[delta - 1],
                      "valid branching coefficient must remain unchanged");
            CPhotoDisintegrationCanal single(isotope, ESingle);
            CPhotoDisintegrationCanal doubly(isotope, EDouble);
            for (int delta = -1; delta <= 4; ++delta) {
                check(single.getRate(delta) == (delta == 1 ? 1. : 0.), "single channel");
                check(doubly.getRate(delta) == (delta == 2 ? 1. : 0.), "double channel");
            }
        }
    } catch (Utils::Exception* error) {
        std::cerr << error->Message() << std::endl;
        delete error;
        return 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}
