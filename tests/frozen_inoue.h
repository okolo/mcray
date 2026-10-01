#pragma once
#include "Inoue12IROSpectrum.h"

// A test fixture for the agreed static-medium experiment; production EBL
// evolution remains supplied by the existing background implementation.
class FrozenInoueBaseline : public Backgrounds::Inoue12BaselineIROSpectrum {
public:
    double n(double energy, double) const override {
        return Backgrounds::Inoue12BaselineIROSpectrum::n(energy, 0.);
    }
};
