/*
 * Deflection3D.cpp
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
#include "Deflection3D.h"
#include "Test.h"


namespace Interactions {
    using namespace mcray;
    using namespace Utils;


    Deflection3D::Deflection3D(MagneticField * aMagnField, u_int aAccuracy) :
            fGlobalB(aMagnField),
            fAccuracy(aAccuracy)
    {
    }

    Deflection3D::Deflection3D(uint aAccuracy, int aUniqueMFslot):
            fAccuracy(aAccuracy),
            fMFslot(aUniqueMFslot)
    {
        if(fMFslot<0)
            fMFslot = Particle::ReserveInteractionDataSlot();
        if(fMFslot<0)
            Exception::Throw("Deflection3D: Failed to reserve interaction slot");
    }

    void Deflection3D::GetRotationRate(const Particle &aParticle, const std::vector<double> &aBgauss, std::vector<double>& outRate)
    {
        //dn/dt = qeB/E, where q is charge in units of electron charge; e=sqrt(alpha)-electron charge
        double mult = aParticle.ElectricCharge()*units.Gauss*units.e/aParticle.Energy;
        const double* N = aParticle.Pdir;
        ASSERT(outRate.size()==3);
        outRate[0] = (N[1]* aBgauss[2]- aBgauss[1]*N[2])*mult;
        outRate[1] = (N[2]* aBgauss[0]- aBgauss[2]*N[0])*mult;
        outRate[2] = (N[0]* aBgauss[1]- aBgauss[0]*N[1])*mult;
    }

    double Deflection3D::Rate(const Particle &aParticle) const {
        if(!aParticle.ElectricCharge())
            return 0;
        std::vector<double> deflRate(3);
        std::vector<double> curB(3);
        MagneticField* B = fGlobalB;
        if(!B)
            B = (MagneticField*)(aParticle.interactionData[fMFslot]);
        ASSERT(B);
        B->GetValueGauss(aParticle.X, aParticle.Time, curB);
        GetRotationRate(aParticle, curB, deflRate);
        double absRate = sqrt(deflRate[0]*deflRate[0]+deflRate[1]*deflRate[1]+deflRate[2]*deflRate[2]);
        return absRate;
    }

    bool Deflection3D::Propagate(cosmo_time deltaT, Particle &aParticle, Randomizer &aRandomizer) const {
        if(!aParticle.ElectricCharge())
            return false;
        cosmo_time initialT = aParticle.Time.t();
        cosmo_time beta = aParticle.beta();//particle speed
        MagneticField* pB = fGlobalB;
        if(!pB)
            pB = (MagneticField*)(aParticle.interactionData[fMFslot]);
        ASSERT(pB);
        cosmo_time dx = pB->MinVariabilityScale(aParticle.Time)/fAccuracy;
        u_int nSteps = (u_int)(deltaT/dx + 1.);
        dx = deltaT/nSteps;
        std::vector<double> deflRate(3);
        std::vector<double> curB(3);
        double newPdir[3];
        int totSteps=0;//debug info

        for(u_int stepB = 0; stepB<nSteps; stepB++){
            pB->GetValueGauss(aParticle.X, aParticle.Time, curB);
            GetRotationRate(aParticle, curB, deflRate);
            double absRate = sqrt(deflRate[0]*deflRate[0]+deflRate[1]*deflRate[1]+deflRate[2]*deflRate[2]);
            cosmo_time maxStep = 0.05*M_PI/absRate/fAccuracy;//max step corresponds to 9/fAccuracy degrees turn
            u_int nSubSteps = (u_int)(dx/maxStep+1.);
            cosmo_time dXsubStep = dx/nSubSteps;
            for(u_int substep=0; substep<nSubSteps; substep++){
                if(substep)
                    GetRotationRate(aParticle, curB, deflRate);
                double absPdir = 0.;
                for(int i=0; i<3; i++){
                    double val = aParticle.Pdir[i] + deflRate[i]*dXsubStep;
                    newPdir[i] = val;
                    absPdir += val*val;
                }
                absPdir = sqrt(absPdir);//norm of new vector
                for(int i=0; i<3; i++){
                    double val = newPdir[i]/absPdir;//make sure rotated vector has unit norm
                    aParticle.X[i] += 0.5*(aParticle.Pdir[i] + val)*beta*dXsubStep;
                    aParticle.Pdir[i] = val;
                }
                aParticle.Time += dXsubStep;
                totSteps++;
            }
        }
        //double timeError = fabs((aParticle.Time.t()-initialT-deltaT)/deltaT)/totSteps;
        //double eps = deltaT/initialT;
        //char a[1024];
        //sprintf(a, "eps=%e,timeError=%e", eps, timeError);
        //ASSERT(timeError<1e-14); # adding small dt to large t leads to calculation error. We minimize the error by using long double for time
        aParticle.Time.setT(initialT + deltaT);
        //aParticle.Time.setT(finalT);//increase accuracy
        return true;
    }

    DeflectionInteraction *Deflection3D::Clone() const {
        return fGlobalB ? (new Deflection3D(fGlobalB, fAccuracy)) : (new Deflection3D(fAccuracy, fMFslot));
    }

    MonochromaticMF::MonochromaticMF(Randomizer &aRandomizer, double aLamdbaMpc,
                                                       double aAmplitudeGauss):
            fLambda(aLamdbaMpc*units.Mpc),
            fK(2.*M_PI/fLambda),
            fAbsBgauss(aAmplitudeGauss),
            fBeta(aRandomizer.Rand()*M_PI*2.),
            fAlpha(aRandomizer.Rand()*M_PI*2.),
            fAmplitude(3)
    {
        /// Choosing random direction: cos(Theta) and phi are distributed uniformly
        fCosTheta = 1.-(aRandomizer.Rand()*2);
        fSinTheta = sqrt(1.- fCosTheta * fCosTheta);
        double phiK(aRandomizer.Rand()*M_PI*2.);
        fCosPhi = cos(phiK);
        fSinPhi = sin(phiK);
        init();
    }

    MonochromaticMF::MonochromaticMF(double aLamdbaMpc, double aAmplitudeGauss, double aBetaK,
    double aAlphaK, double aTheta, double aPhi):
            fLambda(aLamdbaMpc*units.Mpc),
            fK(2.*M_PI/fLambda),
            fAbsBgauss(aAmplitudeGauss),
            fBeta(aBetaK),
            fAlpha(aAlphaK),
            fAmplitude(3),
            fCosTheta(cos(aTheta)),
            fSinTheta(sin(aTheta)),
            fCosPhi(cos(aPhi)),
            fSinPhi(sin(aPhi))
    {
        init();
    }

    void MonochromaticMF::init()
    {
        double cosAlpha = cos(fAlpha);
        double sinAlpha = sin(fAlpha);

        fAmplitude[0] = std::complex<double>(cosAlpha* fCosTheta * fCosPhi, -sinAlpha* fSinPhi)* fAbsBgauss;
        fAmplitude[1] = std::complex<double>(cosAlpha* fCosTheta * fSinPhi, sinAlpha* fCosPhi)* fAbsBgauss;
        fAmplitude[2] = std::complex<double>(-cosAlpha* fSinTheta, 0)* fAbsBgauss;
    }

    void MonochromaticMF::GetValueGauss(const double *x, const CosmoTime &aTime,
                                               std::vector<double> &outValue) const
    {
        ASSERT(outValue.size()==3);
        double zPrime = fSinTheta * fCosPhi *x[0] + fSinTheta * fSinPhi *x[1] + fCosTheta *x[2];
        std::complex<double> phase(0,fK*zPrime+ fBeta);
        std::complex<double> mult = exp(phase);
        for(int i=0; i<3; i++){
            outValue[i] = (fAmplitude[i]*mult).real();
        }
    }

    double MonochromaticMF::MinVariabilityScale(const CosmoTime &aTime) const {
        return fLambda;
    }

    int MagneticField::Print(double lambdaMpc, std::ostream& aOutput) {
        std::vector<double> B(3);
        CosmoTime t(0);
        double step=lambdaMpc*units.Mpc/20.;
        double L=lambdaMpc*units.Mpc*3;
        for(double x=0.; x<=L; x+=step) {
            for (double y = 0.; y <= L; y+=step) {
                for (double z = 0.; z <= L; z += step) {
                    const double r[] = {x, y, z};
                    GetValueGauss(r, t, B);
                    aOutput << x / units.Mpc << "\t" << y / units.Mpc << "\t" << z / units.Mpc << "\t"
                              << B[0] << "\t" << B[1] << "\t" << B[2] << "\n";
                }
                aOutput << "\n";
            }
            aOutput << "\n\n";
        }
        aOutput << "# use gnuplot command \"set pm3d at b; splot 'B' i 0 u 2:3:col\"  where col=4,5 or 6 to view the data\n";
        return 0;
    }

    int MonochromaticMF::UnitTest() {
        std::vector<double> B(3);
        CosmoTime t(0);
        Randomizer r;
        double lambdaMpc = 1.;
        double Bgauss = 10.;
        MonochromaticMF mf(r,lambdaMpc,Bgauss);
        mf.Print(lambdaMpc);
    }

    int Deflection3D::UnitTest() {
        if(!cosmology.IsInitialized())
            cosmology.Init();
        CosmoTime startTime(0.01);
        Particle e(Electron, startTime.z());
        e.Energy = 1e13*units.eV;
        Randomizer r;
        double lambdaMpc = 1;
        double Bgauss = 1e-15;
        double dT = units.Mpc*1.;

        double Beta = 0.1;
        double Alpha = M_PI/4;
        double Theta = 0;
        double Phi = 0;
        SmartPtr<MagneticField> mf = new MonochromaticMF(lambdaMpc, Bgauss, Beta, Alpha, Theta, Phi);
        std::vector<double> B0(3);
        mf->GetValueGauss(e.X, e.Time, B0);
        double Bperp = sqrt(B0[0]*B0[0]+B0[1]*B0[1]);
        Deflection3D defl(mf);
        //Deflection3D defl(new MonochromaticMF(r,lambdaMpc,Bgauss));
        defl.Propagate(dT, e, r);
        double L = sqrt(e.X[0]*e.X[0]+e.X[1]*e.X[1]+e.X[2]*e.X[2]);
        double P = sqrt(e.Pdir[0]*e.Pdir[0]+e.Pdir[1]*e.Pdir[1]+e.Pdir[2]*e.Pdir[2]);
        double theta = acos(e.Pdir[2])/M_PI*180.;

        std::cout << e.X[0]/units.Mpc << "\t" << e.X[1]/units.Mpc << "\t" << e.X[2]/units.Mpc << "\n";
        std::cout << "B(0)=" << Bperp << "\t(dt-dL)/dt = " << (dT-L)/dT << "\t|Pdir|=" << P << "\ttheta=" << theta << " deg.\n";
        return 0;
    }

    void TurbulentMF::GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const {
        std::vector<double> wave(3, 0.);
        outValue.assign(3, 0.);
        for(std::vector<SafePtr<MonochromaticMF> >::const_iterator pw = fWaves.begin(); pw!=fWaves.end(); pw++){
            (*pw)->GetValueGauss(x, aTime, wave);
            for(uint i=0; i<3; i++)
                outValue[i] += std::sqrt(2)*(wave[i]);
        }
    }

    double TurbulentMF::MinVariabilityScale(const CosmoTime &aTime) const {
        return fLmax*0.02;
    }

    TurbulentMF::TurbulentMF(Randomizer &aRandomizer, double aLmin_Mpc, double aLmax_Mpc, double aVariance_Gauss, int aNmodes):
            fLmin(aLmin_Mpc*units.Mpc),
            fLmax(aLmax_Mpc*units.Mpc)
    {
        ASSERT(aNmodes>1.);
        ASSERT(aLmax_Mpc>aLmin_Mpc);
        const double gamma = 11./3.;
        std::vector<double> norms;
        std::vector<double> Ks;
        double sumNorm = 0;
        double aMultK = pow(aLmax_Mpc/aLmin_Mpc, 1./(aNmodes-1.));
        double k = 2*M_PI/fLmax;
        // std::cout << aLmin_Mpc << '\t' << aLmax_Mpc << '\t' << aNmodes << '\t' << aMultK << '\n';
        for(uint i=0; i<aNmodes; i++){
            Ks.push_back(k);
            double normK = (k*fLmax)*(k*fLmax)*(k*fLmax)*pow(k*fLmax,-gamma);
            norms.push_back(normK);
            sumNorm += normK;
            fWaves.push_back(0);
            k *= aMultK;
        }
        for(uint i=0; i<norms.size(); i++){
            fWaves[i] = new MonochromaticMF(aRandomizer, 2.*M_PI/Ks[i]/units.Mpc, aVariance_Gauss*sqrt(norms[i]/sumNorm));
        }
    }

    int TurbulentMF::UnitTest() {
        if(!cosmology.IsInitialized())
            cosmology.Init();
        CosmoTime startTime(0.03);
        Particle e(Electron, startTime.z());
        e.Energy = 1e13*units.eV;
        Randomizer r;
        double lambdaMpc = 1;
        double Bgauss = 1e-15;
        SmartPtr<MagneticField> mf = new TurbulentMF(r, lambdaMpc/20, lambdaMpc*5, Bgauss, 100);
        Deflection3D defl(mf);
        std::ofstream outB("TurbulentMF_B.txt");
        mf->Print(lambdaMpc, outB);
        std::vector<double> B0(3);
        double t=0;
        std::cout << "E/eV = " << e.Energy/units.eV << "\tB/G = " << Bgauss << "\tLc/Mpc = " << lambdaMpc << "\n";
        mf->GetValueGauss(e.X, e.Time, B0);
        double Bperp = sqrt(B0[0]*B0[0]+B0[1]*B0[1]);
        std::cout << "t/Mpc=0\tB_xy/G=" << Bperp << "\tr = (" << e.X[0]/units.Mpc << "\t" << e.X[1]/units.Mpc << "\t" << e.X[2]/units.Mpc << ")\n";

        for(double tMpc=1.; tMpc<100; tMpc*=3) {
            defl.Propagate(units.Mpc * (tMpc - t), e, r);
            mf->GetValueGauss(e.X, e.Time, B0);
            double Bperp = sqrt(B0[0]*B0[0]+B0[1]*B0[1]);
            std::cout << "t/Mpc=" << tMpc << "\tB_xy/G=" << Bperp << "\tr = (" << e.X[0]/units.Mpc << "\t" << e.X[1]/units.Mpc << "\t" << e.X[2]/units.Mpc << ")\n";
            double L = sqrt(e.X[0]*e.X[0]+e.X[1]*e.X[1]+e.X[2]*e.X[2]);
            double P = sqrt(e.Pdir[0]*e.Pdir[0]+e.Pdir[1]*e.Pdir[1]+e.Pdir[2]*e.Pdir[2]);
            double theta = acos(e.Pdir[2])/M_PI*180.;
            std::cout << "(dt-dL)/dt = " << (tMpc-L/units.Mpc)/tMpc << "\t|Pdir| = " << P << "\ttheta = " << theta << " deg.\n";
            t=tMpc;
        }

        return 0;
    }
}