#pragma once
#include "observer_api.h"
#include "transport_api.h"
struct CoupledInput {
    ObserverInput observer;
    int particle; // electron=0, positron=1, proton=2, photon=3
    int accuracy,fieldKind; // free, constant, saved grid
    double energy;
};
struct CoupledOutput {
    TransportOutput transport;
    bool hit,pastEnd;
    long double endError;
    unsigned long moves;
};
CoupledOutput HistoricalCoupled(const CoupledInput&,const char* grid,const double* field,double scale);

// An immutable field may be loaded once per independent source stream.
void* HistoricalCoupledField(int kind,const char* grid,const double* field,double scale);
void DeleteHistoricalCoupledField(void*);
CoupledOutput HistoricalCoupledWithField(const CoupledInput&,void* field);
