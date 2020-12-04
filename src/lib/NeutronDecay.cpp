/*
 * NeutronDecay.cpp
 *
 * Author:
 *       Oleg Kalashev
 *
 * Copyright (c) 2020 Institute for Nuclear Research, RAS
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


#include "NeutronDecay.h"
#include "Test.h"

namespace Interactions {
    using namespace mcray;

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
