/*
 * PhotoDisintegration.h
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


#ifndef MCRAY_PHOTODISINTEGRATION_H
#define MCRAY_PHOTODISINTEGRATION_H

#include "Interaction.h"
#include "Background.h"
#include "Nucleus.h"

namespace Interactions {
    using namespace mcray;

    class PhotoDisintegration : public RandomInteraction {
        class SigmaS : public Function
        {
        public:
            SigmaS(const CPDMLine* canals);
            double f(double s) const;
            double Xmin() const {return fXmin;}
            double Xmax() const {return fXmax;}
        private:
            const CPDMLine* fCanals;
            double fXmin;
            double fXmax;
            double fM;
        };
        class SigmaK : public Function
        {
        public:
            SigmaK(const CPDMLine* canals);
            double f(double s) const;
            double Xmin() const {return fXmin;}
            double Xmax() const {return fXmax;}
        private:
            const CPDMLine* fCanals;
            double fXmin;
            double fXmax;
        };
    public:
        PhotoDisintegration(BackgroundIntegral* aBackground, bool aSampleK=true);
        RandomInteraction* Clone() const;
        virtual ~PhotoDisintegration();
        double Rate(const Particle& aParticle) const;
        static void UnitTest(ParticleType aPrimary, bool aPrintRates, bool aPropagate, bool aSplit, double aEmin, double aEmax, double aZmax, int nParticles=1000);
        static int UnitTest();
        bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
        bool GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const;
        void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aL, Randomizer& aRandomizer) const;
        bool SampleK(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
        void SampleSecondariesWithS(Particle &aParticle, std::vector<Particle> &aSecondaries, double aS, Randomizer &aRandomizer) const;
        void SampleSecondariesWithK(Particle &aParticle, std::vector<Particle> &aSecondaries, double aK, Randomizer &aRandomizer) const;

    private:
        CPhotoDisintegrationMap         fMap;
        SmartPtr<BackgroundIntegral>	fBackground;
        bool                            fSampleK;
    };
}

#endif //MCRAY_PHOTODISINTEGRATION_H
