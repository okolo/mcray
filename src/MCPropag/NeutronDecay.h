//
// Created by ok on 01/08/2015.
//

#ifndef MCPROPAG_NEUTRONDECAY_H
#define MCPROPAG_NEUTRONDECAY_H

#include "Interaction.h"

namespace Interactions {
    using namespace MCPropag;

    class NeutronDecay : public RandomInteraction {

    public:
        NeutronDecay();

        virtual double Rate(const Particle &aParticle) const;

        virtual bool GetSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries,
                                    Randomizer &aRandomizer) const;

        virtual RandomInteraction* Clone() const;
        static int UnitTestEconserv();
    private:
        static const double Neutron_lifetime; // neutron lifetime in it's rest frame (sec)
        const double n_decayM; // neutron decay probability it's rest frame times neutron mass
        double n_decay_k; // secondary leptons momenta in CM frame (assuming zero nu mass)
        double n_decay_E; // electron energy in CM frame
    };
}
#endif //MCPROPAG_NEUTRONDECAY_H
