/*
 * Deflection3D.h
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

#ifndef MCRAY_DEFLECTION3D_H
#define MCRAY_DEFLECTION3D_H

#include <complex>
#include "Interaction.h"
#include "MathUtils.h"


namespace Interactions {
    using namespace mcray;
    using namespace Utils;

    class MagneticField : public SmartReferencedObj{
    public:
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double>& outValue) const = 0;

        virtual double MinVariabilityScale(const CosmoTime &aTime) const = 0;

        int Print(double lambdaMpc, std::ostream& aOutput=std::cout);

        double CorrelationLength(std::vector<double> x, std::vector<double> aDirection, const CosmoTime &aTime,
                                         double aMaxValMpc=10000);

    };

    class Deflection3D : public DeflectionInteraction {
    public:
        //aMagnField should provide magnetic field in Gauss
        Deflection3D(MagneticField * aMagnFieldGauss, uint aAccuracy=10);

        //use custom magnetic field from each particle
        //@param aUniqueMFslot custom data slot to use or -1 to reserve one
        Deflection3D(uint aAccuracy = 10, int aUniqueMFslot = -1);

        virtual double Rate(const Particle &aParticle) const;

        virtual bool Propagate(cosmo_time deltaT, Particle &aParticle, Randomizer &aRandomizer) const;

        virtual DeflectionInteraction *Clone() const;

        //return 3D vector dn/dt === 1/beta dV/dt
        static void GetRotationRate(const Particle &aParticle, const std::vector<double> &aBgauss, std::vector<double>& outRate);
        static int UnitTest();
        int MFSlot() const { return fMFslot; }
    private:
        SmartPtr<MagneticField> fGlobalB;
        u_int fAccuracy;
        int fMFslot;
    };

    //Randomly oriented plane wave with given k
    //Generated as fixed k term in formula (3) of J. Giacalone, J.R. Jokipii, 1999, Ap.J., 520:204-214
    class MonochromaticMF : public MagneticField{
    public:
        MonochromaticMF(Randomizer &aRandomizer, double aLamdbaMpc, double aAmplitudeGauss);
        //constructor with specific orientation and phase
        MonochromaticMF(double aLamdbaMpc, double aAmplitudeGauss, double aBetaK,
                        double aAlphaK, double aTheta, double aPhi);

        virtual void GetValueGauss(const double *x, const CosmoTime &aTime,
                                   std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;

        static int UnitTest();
    private:
        void init();
    private:
        double fLambda;
        double fK;//2pi/lambda
        double fAbsBgauss;
        double fBeta;
        double fAlpha;
        double fCosTheta;
        double fSinTheta;
        double fCosPhi;
        double fSinPhi;

        std::vector<std::complex<double> >fAmplitude;
    };

    //Turbulent magnetic field defined as sum of randomly oriented plane waves
    //as in J. Giacalone, J.R. Jokipii, 1999, Ap.J., 520:204-214 formulas (3)-(7)
    class TurbulentMF : public MagneticField{

    public:
        TurbulentMF(Randomizer &aRandomizer, double aLmin_Mpc, double aLmax_Mpc, double aVariance_Gauss, int aNmodes);
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
        static int UnitTest();
        static double MeanCorLength(Randomizer &aRandomizer, double aLmin_Mpc, double aLmax_Mpc, double aVariance_Gauss,
                                    int aNmodes, const CosmoTime &aTime, int aNsamples=1000, double aMaxValMpc=10000);
    private:
        double fLmin;
        double fLmax;
        std::vector<SafePtr<MonochromaticMF> > fWaves;
    };

}
#endif //MCRAY_DEFLECTION3D_H
