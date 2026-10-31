//
// Created by ok on 25/03/2016.
//

#ifndef MCPROPAG_PHOTODISINTEGRATION_H
#define MCPROPAG_PHOTODISINTEGRATION_H

#include "Interaction.h"
#include "Background.h"
#include "Nucleus.h"

namespace Interactions {
    using namespace MCPropag;

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

#endif //MCPROPAG_PHOTODISINTEGRATION_H
