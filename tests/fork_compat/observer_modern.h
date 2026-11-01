#pragma once
#include "observer_api.h"
#include <memory>
namespace mcray { class Particle; class SphereObserver; }
namespace observer_validation {
std::shared_ptr<mcray::SphereObserver> Observer(const ObserverInput&);
mcray::Particle Initial(const ObserverInput&,double weight=2.5);
ObserverFlight ModernFlight(const ObserverInput&,bool variable,double weight=2.5);
void PrintFlight(const ObserverFlight&);
}
