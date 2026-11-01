#pragma once
#include "observer_api.h"
#include "transport_api.h"
struct CoupledInput {
    ObserverInput observer;
    int particle,accuracy,fieldKind; // free, constant, saved grid
    double energy;
};
struct CoupledOutput {
    TransportOutput transport;
    bool hit,pastEnd;
    long double endError;
    unsigned long moves;
};
CoupledOutput HistoricalCoupled(const CoupledInput&,const char* grid,const double* field,double scale);
