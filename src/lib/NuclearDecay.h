/*
 * NuclearDecay.h
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


#ifndef MCRAY_NUCLEARDECAY_H
#define MCRAY_NUCLEARDECAY_H

#include "Interaction.h"

namespace Interactions {
    using namespace mcray;

    bool UnstableDecay(int aPrimZ, int aPrimN, Particle &aParticle, std::vector<Particle> &aSecondaries,
                            Randomizer &aRandomizer);

    inline std::string                     TableLifetime = "tables/nuclear_decay.txt";
    inline int                             MaxNBranches = 32;

    class NuclearDecay : public RandomInteraction {

    public:
        NuclearDecay();

        virtual double Rate(const Particle &aParticle) const;

        virtual bool GetSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries,
                                    Randomizer &aRandomizer) const;
                 
        virtual RandomInteraction* Clone() const;

    };
}
#endif 
