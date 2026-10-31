//
// Created by ok on 01/08/2015.
//

#include "NeutronDecay.h"
#include "Test.h"

namespace Interactions {
    using namespace MCPropag;

    const double NeutronDecay::Neutron_lifetime=881.5; // (sec)
    NeutronDecay::NeutronDecay():
            n_decayM(Particle::Mass(Neutron)/Neutron_lifetime/units.second)
    {
        double Me = Particle::Mass(Electron);
        double delta_m = Particle::Mass(Neutron)-Particle::Mass(Proton);
        n_decay_k=0.5*(delta_m-Me*Me/delta_m);
        n_decay_E=sqrt(n_decay_k*n_decay_k+Me*Me);
    }

    double NeutronDecay::Rate(const Particle &aParticle) const {
        if(aParticle.Type!=Neutron)
            return 0;
        return n_decayM/aParticle.Energy;
    }

    bool NeutronDecay::GetSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries,
                                      Randomizer &aRandomizer) const {
        if(aParticle.Type!=Neutron)
            return false;

        //proton
        Particle proton = aParticle;
        proton.Type = Proton;
        proton.Energy *= Particle::Mass(Proton)/Particle::Mass(Neutron);
        double gamma = aParticle.Energy/aParticle.Mass();
        double beta = aParticle.beta();

        //electron
        //sample electron momentum direction in CM frame
        double angle_e = aRandomizer.Rand()*M_PI;
        double cos_e = cos(angle_e);
        Particle electron = aParticle;
        electron.Type = Electron;
        electron.Energy = gamma*(n_decay_E-beta*cos_e*n_decay_k);

        //neutrino
        Particle nu = aParticle;
        nu.Type = NeutrinoAE;
        nu.Energy = gamma*n_decay_k*(1.0+beta*cos_e);//cos_nu=-cos_e

        aSecondaries.push_back(proton);
        aSecondaries.push_back(electron);
        aSecondaries.push_back(nu);

        return true;
    }

    int NeutronDecay::UnitTestEconserv(){
        double aZ=0.1;
        double aE=1e19;
        int nParticles = 100;
        Test::EConservationTest test(aZ,1e6,2015);
        test.alphaThinning = 0.0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
        //test.Engine().AddInteraction(new GZK(test.Backgr(),(int)(test.GetRandomizer().CreateIndependent())));
        test.Engine().AddInteraction(new NeutronDecay());
        Particle p(Neutron, aZ);
        p.Energy = aE*units.eV;
        test.SetPrimary(p, nParticles);
        test.Run();
        return 0;
    }

    RandomInteraction *NeutronDecay::Clone() const {
        return new NeutronDecay();
    }
}
