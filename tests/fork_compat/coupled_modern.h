#pragma once
#include "coupled_api.h"
#include "Deflection3D.h"

// The same modern path serves fixed-state controls and the source cost pilot.
Utils::SmartPtr<Interactions::MagneticField> ModernCoupledField(
    int kind,const char* grid,const double* field,double scale);
CoupledOutput ModernCoupled(const CoupledInput&,Interactions::MagneticField&,double weight=1);
void PrintCoupled(const CoupledOutput&);
