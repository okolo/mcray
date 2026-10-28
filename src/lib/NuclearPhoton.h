#ifndef MCRAY_NUCLEAR_PHOTON_H
#define MCRAY_NUCLEAR_PHOTON_H
#include "Interaction.h"
namespace Interactions {
struct NuclearPhotonLine { double energy=0., probability=0., rawYield=0.; };
inline void EmitNuclearPhoton(const mcray::Particle& emitter,double restEnergy,
                              std::vector<mcray::Particle>& products,mcray::Randomizer& random) {
    // Approved CRPropa ultrarelativistic boost; nuclear energy is not debited.
    mcray::Particle photon=emitter;
    photon.SetType(mcray::Photon);
    photon.Energy=restEnergy*(emitter.Energy/emitter.Mass())*(1.-(2.*random.Rand()-1.));
    if (!std::isfinite(photon.Energy) || photon.Energy<0.) Utils::Exception::ThrowFatal("Invalid boosted nuclear photon energy");
    photon.fCascadeProductionTime=emitter.Time;
    products.push_back(photon);
}
}
#endif
