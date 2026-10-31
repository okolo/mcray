//
// Created by ok on 06/08/2015.
//

#ifndef MCPROPAG_DEFLECTION3D_H
#define MCPROPAG_DEFLECTION3D_H

#include <complex>
#include <string>
#include "Interaction.h"
#include "MathUtils.h"


namespace Interactions {
    using namespace MCPropag;
    using namespace Utils;

    class MagneticField : public SmartReferencedObj{
    public:
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double>& outValue) const = 0;

        virtual double MinVariabilityScale(const CosmoTime &aTime) const = 0;

        int Print(double lambdaMpc, std::ostream& aOutput=std::cout);
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
        TurbulentMF(Randomizer &aRandomizer, double aLcor_Mpc, double aVariance_Gauss, 
                    double aMultK=1.05/* steps in K */, double aMinVarScale=0.1);
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
        static int UnitTest();
    private:
        double fLc;
        double fMinVarScale;
        std::vector<SafePtr<MonochromaticMF> > fWaves;
    };


    class CellMF : public MagneticField{
    public:
        CellMF(Randomizer &aRandomizer, double aLcor_Mpc, double aAmplitudeGauss);
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
    private:
        double fLc;
        double fAmplitude;
        double fCurrentMaxz;
        int fSize;
        std::vector<int> fCellDirection;
    };

    class TableLOS_MF : public MagneticField{
    public:
        TableLOS_MF(std::string filename);
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
    private:
        std::vector<double> fL;
        std::vector<double> fBx, fBy, fBz;
        int fnL;
    };


    class TableLOSTube_MF : public MagneticField{
    public:
        TableLOSTube_MF();
        virtual void GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
    private:
        std::vector<double> fL;
        std::vector< std::vector<double> > fBx, fBy, fBz;
        int fnL;
    };

}
#endif //MCPROPAG_DEFLECTION3D_H
