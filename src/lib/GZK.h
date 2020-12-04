/*
 * GZK.h
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

#ifndef MCRAY_GZK_H
#define MCRAY_GZK_H

#include "Interaction.h"
#include "Background.h"

namespace Interactions {
    using namespace mcray;

    class GZK : public RandomInteractionS {
    public:
        GZK(BackgroundIntegral* aBackground, int aRandSeed=0);

        virtual double Rate(const Particle &aParticle) const;

        virtual RandomInteraction *Clone() const;

        virtual bool SampleS(const Particle &aParticle, double &aS, Randomizer &aRandomizer) const;

        virtual void SampleSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries, double aS,
                                       Randomizer &aRandomizer) const;
        static void UnitTest();
        static void UnitTestEconserv(double aE, double aZ);
    private:
        Function* InitSigma(ParticleType aPrim);
        SmartPtr<BackgroundIntegral>	fBackground;
        SafePtr<Function>               fSigmaN;
        SafePtr<Function>               fSigmaP;
        double                          fMpSophia;
        double                          fMnSophia;
    };
}

#endif //MCRAY_GZK_H
