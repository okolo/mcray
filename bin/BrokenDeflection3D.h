//
// Created by ok on 06/08/2015.
//

#ifndef MCPROPAG_DEFLECTION3D_H
#define MCPROPAG_DEFLECTION3D_H

#include <complex>
#include <string>
#include <cmath>
#include "Interaction.h"
#include "MathUtils.h"

extern "C" {
#include <fftw3.h>
}

namespace Interactions {
    using namespace MCPropag;
    using namespace Utils;

    class MagneticField : public SmartReferencedObj{
    public:
        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double>& outValue) const = 0;

        virtual double MinVariabilityScale(const CosmoTime &aTime) const = 0;

        int Print(double stepMpc,double maxMpc, std::ostream& aOutput=std::cout);
        double get_div(const double x_Mpc, const double y_Mpc, const double z_Mpc, const double delta_Mpc=0.001);
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

    // Load field from file
    class SavedMF : public MagneticField{
    public:
        SavedMF(const char* fname);
        SavedMF(){Nx=0;Ny=0;Nz=0;}; // This constructor must be overloaded
        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime,
                                   std::vector<double> &outValue) const;
        virtual void GetValueGauss_NEAREST(const long double *x, const CosmoTime &aTime,
                                   std::vector<double> &outValue) const;
        virtual void GetValueGauss_LINEAR(const long double *x, const CosmoTime &aTime,
                                   std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
        ~SavedMF();
    private:
    	size_t coord2indexZero(double coord, double min_val, double max_val, size_t N) const;
    	size_t coord2index(double x, double y, double z) const;
    	double linApprox(const double* f,
    					size_t i, size_t j, size_t k,
    					double x, double y, double z) const;
    protected:
        double* B_X;
        double* B_Y;
        double* B_Z;
        double x_min;
        double y_min;
        double z_min;
        double x_max;
        double y_max;
        double z_max;
        size_t Nx;
        size_t Ny;
        size_t Nz;
    };

    //Randomly oriented plane wave with given k
    //Generated as fixed k term in formula (3) of J. Giacalone, J.R. Jokipii, 1999, Ap.J., 520:204-214
    class MonochromaticMF : public MagneticField{
    public:
        MonochromaticMF(Randomizer &aRandomizer, double aLamdbaMpc, double aAmplitudeGauss);
        //constructor with specific orientation and phase
        MonochromaticMF(double aLamdbaMpc, double aAmplitudeGauss, double aBetaK,
                        double aAlphaK, double aTheta, double aPhi);

        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime,
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
                    double aMultK=1.05/* steps in K */, double aMinVarScale=0.1, double PMFCut_p=-1.);
        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
        static int UnitTest();
    private:
        double fLc;
        double fMinVarScale;
        double PMFCutDistance;
        std::vector<SafePtr<MonochromaticMF> > fWaves;
    };

    //Turbulent magnetic field defined as reversed Fourier from k space
    //as in CRpropa
    class FFT_TurbulentMF : public SavedMF{
    public:
        FFT_TurbulentMF(Randomizer &aRandomizer, double aLcor_Mpc, double aVariance_Gauss,
                    double cube_size, double aMinVarScale=0.1);
//        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
//        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
//    private:
//        double fLc;
//        double fMinVarScale;
//        double PMFCutDistance;
//        std::vector<SafePtr<MonochromaticMF> > fWaves;
    };


    class CellMF : public MagneticField{
    public:
        CellMF(Randomizer &aRandomizer, double aLcor_Mpc, double aAmplitudeGauss);
        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
    private:
        double fLc;
        double fAmplitude;
        double fCurrentMaxz;
        int fSize;
        std::vector<int> fCellDirection;
    };


    //Galactic magnetic field
    //See: Pshirkov, Tinyakov, Kronberg Newton-McGee 2011 - DOI: 10.1088/0004-637X/738/2/192, arXiv:1103.0814
    class GalacticMF : public MagneticField{
    public:
        GalacticMF();
        virtual void GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const;
        virtual double MinVariabilityScale(const CosmoTime &aTime) const;
        static int UnitTest();
    private:
		double pitch, cos_pitch, sin_pitch, PHI, cos_PHI;  // pitch angle parameters
		double d;     // distance to first field reversal
		double R_sun; // distance between sun and galactic center
		double R_c;   // radius of central region
		double z0_D;    // vertical thickness in the galactic disk
		double B0_D;    // magnetic field scale

		// halo parameters
		double z0_H;  // halo vertical position
		double R0_H;  // halo radial position
		double B0_Hn; // halo magnetic field scale (north)
		double B0_Hs; // halo magnetic field scale (south)
		double z11_H; // halo vertical thickness towards disc
		double z12_H; // halo vertical thickness off the disk
    };
}
#endif //MCPROPAG_DEFLECTION3D_H
