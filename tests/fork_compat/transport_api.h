#pragma once

// Shared POD boundary only; historical types stay in a separate translation unit.
struct TransportInput {
    int particle, accuracy, outerSteps;
    double energy;
    long double position[3], direction[3], initialTime, duration;
};
struct TransportOutput {
    long double position[3], direction[3], elapsed, beta;
    double mass, charge, energy, initialRate, initialScale, globalScale;
    unsigned long fieldCalls, scaleCalls, zeroScaleCalls;
};
TransportOutput HistoricalTransport(const TransportInput& input, const char* grid,
                                    int mode, const double* field, double scale);
