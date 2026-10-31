//
// Created by ok on 06/08/2015.
//

//#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
//#include <boost/numeric/ublas/storage.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/geometry/geometry.hpp>
#include "Deflection3D.h"
#include "Test.h"

namespace Interactions {
    using namespace MCPropag;
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

    double MagneticField::get_div(const double x_Mpc, const double y_Mpc, const double z_Mpc, const double delta_Mpc) const {
        double B_der_x, B_der_y, B_der_z;
        long double x, y, z;
        x = x_Mpc*units.Mpc;
        y = y_Mpc*units.Mpc;
        z = z_Mpc*units.Mpc;

        long double dx, dy, dz;
        dx = delta_Mpc*units.Mpc;
        dy = delta_Mpc*units.Mpc;
        dz = delta_Mpc*units.Mpc;

        std::vector<double> B1(3), B2(3);
        CosmoTime t(0);
        long double rx1[] = {x + dx, y, z};
        long double rx2[] = {x - dx, y, z};

        long double ry1[] = {x, y + dy, z};
        long double ry2[] = {x, y - dy, z};

        long double rz1[] = {x, y, z + dz};
        long double rz2[] = {x, y, z - dz};

        GetValueGauss(rx1, t, B1);
        GetValueGauss(rx2, t, B2);
        B_der_x = (B1[0] - B2[0])/(2*dx/units.Mpc);

        GetValueGauss(ry1, t, B1);
        GetValueGauss(ry2, t, B2);
        B_der_y = (B1[1] - B2[1])/(2*dy/units.Mpc);

        GetValueGauss(rz1, t, B1);
        GetValueGauss(rz2, t, B2);
        B_der_z = (B1[2] - B2[2])/(2*dz/units.Mpc);
        double divB = B_der_x + B_der_y + B_der_z;

        std::cout << "x = " << x_Mpc << '\n';
        std::cout << "y = " << y_Mpc << '\n';
        std::cout << "z = " << z_Mpc << '\n';
        std::cout << "dBx/dx = " << B_der_x*1e9 << '\t' << "nG/Mpc" << '\n';
        std::cout << "dBy/dy = " << B_der_y*1e9 << '\t' << "nG/Mpc" << '\n';
        std::cout << "dBz/dz = " << B_der_z*1e9 << '\t' << "nG/Mpc" << '\n';
        std::cout << "div B = " << divB*1e9 << '\t' << "nG/Mpc" << '\n';
        std::cout << "-------------" << '\n';
        return divB;
    }

    double MagneticField::get_integral_div(const double x0_Mpc, const double y0_Mpc, const double z0_Mpc, const double x1_Mpc, const double y1_Mpc, const double z1_Mpc, const double step) const {
    	long double x0, y0, z0;
    	long double x1, y1, z1;
    	long double stp = step*units.Mpc;
    	std::cout<<"Integrate with step "<<step<<" Mpc..."<<std::endl;
        x0 = x0_Mpc*units.Mpc;
        y0 = y0_Mpc*units.Mpc;
        z0 = z0_Mpc*units.Mpc;
        x1 = x1_Mpc*units.Mpc;
        y1 = y1_Mpc*units.Mpc;
        z1 = z1_Mpc*units.Mpc;
        long double integralX = 0;
        long double integralY = 0;
        long double integralZ = 0;
    	std::vector<double> B(3);
        CosmoTime t(0);
        //X plates
        for(long double dy=y0;dy<y1;dy+=stp){
	        for(long double dz=z0;dz<z1;dz+=stp){
	        	long double r[] = {x1, dy, dz};
	        	GetValueGauss(r, t, B);
	        	integralX += B[0];
	        	r[0] = x0;
	        	GetValueGauss(r, t, B);
	        	integralX -= B[0];
	        }
        }
        integralX *= step*step;
        std::cout<< "X plates: \t" << integralX*1e9 << " nG*Mpc^2"<< std::endl;
        //Y plates
        for(long double dx=x0;dx<x1;dx+=stp){
	        for(long double dz=z0;dz<z1;dz+=stp){
	        	long double r[] = {dx, y1, dz};
	        	GetValueGauss(r, t, B);
	        	integralY += B[1];
	        	r[1] = y0;
	        	GetValueGauss(r, t, B);
	        	integralY -= B[1];
	        }
        }
        integralY *= step*step;
        std::cout<< "Y plates: \t" << integralY*1e9 << " nG*Mpc^2"<< std::endl;
        //Z plates
        for(long double dx=x0;dx<x1;dx+=stp){
	        for(long double dy=y0;dy<y1;dy+=stp){
	        	long double r[] = {dx, dy, z1};
	        	GetValueGauss(r, t, B);
	        	integralZ += B[2];
	        	r[2] = z0;
	        	GetValueGauss(r, t, B);
	        	integralZ -= B[2];
	        }
        }
        integralZ *= step*step;
        std::cout<< "Z plates: \t" << integralZ*1e9 << " nG*Mpc^2"<< std::endl;
        std::cout<< "Sum of flows:\t" << (integralX+integralY+integralZ)*1e9 << " nG*Mpc^2"<< std::endl;
        std::cout<< "-------------" << std::endl;
        return (integralX+integralY+integralZ);
    }

    void Deflection3D::GetRotationRate(const Particle &aParticle, const std::vector<double> &aBgauss, std::vector<double>& outRate)
    {
        //dn/dt = qeB/E, where q is charge in units of electron charge; e=sqrt(alpha)-electron charge
        double mult = aParticle.ElectricCharge()*units.Gauss*units.e/aParticle.Energy;
        const long double* N = aParticle.Pdir;
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
        cosmo_time dx = pB->MinVariabilityScale(aParticle.X,aParticle.Time,deltaT)/fAccuracy;
        if(dx<=0)dx=deltaT; // zero MinVariabilityScale is stub and must be interpreted as free propagation
        u_int nSteps = (u_int)(deltaT/dx + 1.);
        dx = deltaT/nSteps;
        std::vector<double> deflRate(3);
        std::vector<double> curB(3);
        long double newPdir[3];
        int totSteps=0;//debug info

        for(u_int stepB = 0; stepB<nSteps; stepB++){
            pB->GetValueGauss(aParticle.X, aParticle.Time, curB);
            GetRotationRate(aParticle, curB, deflRate);
            long double absRate = sqrt(deflRate[0]*deflRate[0]+deflRate[1]*deflRate[1]+deflRate[2]*deflRate[2]);
            cosmo_time maxStep = 0.05*M_PI/absRate/fAccuracy;//max step corresponds to 9/fAccuracy degrees turn
            u_int nSubSteps = (u_int)(dx/maxStep+1.);
            cosmo_time dXsubStep = dx/nSubSteps;
            for(u_int substep=0; substep<nSubSteps; substep++){
                if(substep)
                    GetRotationRate(aParticle, curB, deflRate);
                long double absPdir = 0.;
                for(int i=0; i<3; i++){
                    long double val = aParticle.Pdir[i] + deflRate[i]*dXsubStep;
                    newPdir[i] = val;
                    absPdir += val*val;
                }
                absPdir = sqrt(absPdir);//norm of new vector
                for(int i=0; i<3; i++){
                    long double val = newPdir[i]/absPdir;//make sure rotated vector has unit norm
                    aParticle.X[i] += 0.5L*(aParticle.Pdir[i] + val)*beta*dXsubStep;
                    aParticle.Pdir[i] = val;
                }
                aParticle.updateTracer();
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

    void MonochromaticMF::GetValueGauss(const long double *x, const CosmoTime &aTime,
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
        return 0.1*fLambda;
    }

    int MagneticField::Print(double stepMpc, double maxMpc, std::ostream& aOutput) {
        std::vector<double> B(3);
        CosmoTime t(0);
        double step=stepMpc*units.Mpc;
        double L=maxMpc*units.Mpc;
        double BxMean = 0;
        double ByMean = 0;
        double BzMean = 0;
        double BxyMean = 0;
        double BxzMean = 0;
        double ByzMean = 0;
        double Brms = 0;
        double Bmodule = 0;
        aOutput << "#x,Mpc" << "\t" << "y,Mpc" << "\t" << "z,Mpc" << "\t" << "Bx" << "\t" << "By" << "\t" << "Bz" << "\n";

        for(long double x=-L; x<=L; x+=step) {
        	std::cout<<"Saved " << double(int( (L+x)/(2*L) * 1000))/10. << "%...\r";
        	std::cout.flush();
            for (long double y = -L; y <= L; y+=step) {
                for (long double z = -L; z <= L; z += step) {
                    const long double r[] = {x, y, z};
                    GetValueGauss(r, t, B);
                    aOutput << x / units.Mpc << "\t" << y / units.Mpc << "\t" << z / units.Mpc << "\t"
                            << B[0] << "\t" << B[1] << "\t" << B[2] << "\n";
                    BxMean += fabs(B[0]);
                    ByMean += fabs(B[1]);
                    BzMean += fabs(B[2]);
                    BxyMean += sqrt(B[0]*B[0] + B[1]*B[1]);
                    BxzMean += sqrt(B[0]*B[0] + B[2]*B[2]);
                    ByzMean += sqrt(B[1]*B[1] + B[2]*B[2]);
                    Brms += B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
                    Bmodule += sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
                }
                // aOutput << "\n";
            }
            // aOutput << "\n\n";
        }
        std::cout<<"Saved 100%   "<<std::endl;
        BxMean /= pow(2*L/step + 1, 3);
        ByMean /= pow(2*L/step + 1, 3);
        BzMean /= pow(2*L/step + 1, 3);
        BxyMean /= pow(2*L/step + 1, 3);
        BxzMean /= pow(2*L/step + 1, 3);
        ByzMean /= pow(2*L/step + 1, 3);
        Brms /= pow(2*L/step + 1, 3);
        Bmodule /= pow(2*L/step + 1, 3);
        Brms = sqrt(Brms);
        std::cout << "Mean |Bx| =  " << BxMean << '\n';
        std::cout << "Mean |By| =  " << ByMean << '\n';
        std::cout << "Mean |Bz| =  " << BzMean << '\n';
        std::cout << "Mean |Bxy| =  " << BxyMean << '\n';
        std::cout << "Mean |Bxz| =  " << BxzMean << '\n';
        std::cout << "Mean |Byz| =  " << ByzMean << '\n';
        std::cout << "      Brms =  " << Brms << '\n';
        std::cout << "  Mean |B| =  " << Bmodule << '\n';
        // aOutput << "# use gnuplot command \"set pm3d at b; splot 'B' i 0 u 2:3:col\"  where col=4,5 or 6 to view the data\n";
        return 0;
    }

    int MonochromaticMF::UnitTest() {
        std::vector<double> B(3);
        CosmoTime t(0);
        Randomizer r;
        double lambdaMpc = 1.;
        double Bgauss = 10.;
        MonochromaticMF mf(r,lambdaMpc,Bgauss);
        mf.Print(lambdaMpc,lambdaMpc*6.);
        return 0;
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

    void TurbulentMF::GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const {
        std::vector<double> wave(3, 0.);
        outValue.assign(3, 0.);
        if(PMFCutDistance > 0)
        	if( PMFCutDistance*PMFCutDistance < x[0]*x[0]+x[1]*x[1]+x[2]*x[2] )
        		return;
        for(std::vector<SafePtr<MonochromaticMF> >::const_iterator pw = fWaves.begin(); pw!=fWaves.end(); pw++){
            (*pw)->GetValueGauss(x, aTime, wave);
            for(uint i=0; i<3; i++)
                outValue[i] += std::sqrt(2)*(wave[i]);
        }
    }

    double TurbulentMF::MinVariabilityScale(const CosmoTime &aTime) const {
        return fLc*fMinVarScale; // fLc -- типа корр.длина
    }

    TurbulentMF::TurbulentMF(Randomizer &aRandomizer, double aLcor_Mpc, double aVariance_Gauss, double aMultK, double aMinVarScale,double PMFCut_p):
            fLc(aLcor_Mpc*units.Mpc),
            fMinVarScale(aMinVarScale)
    {
        ASSERT(aMultK>1.);
        PMFCutDistance = PMFCut_p*aLcor_Mpc*units.Mpc; // Расстояние обрезания
        uniqHash = "TurbulentMF(";
        uniqHash += std::to_string(aLcor_Mpc)+",";
        uniqHash += std::to_string(aVariance_Gauss)+",";
        uniqHash += std::to_string(aMultK)+",";
        uniqHash += std::to_string(aMinVarScale)+",";
        uniqHash += std::to_string(PMFCut_p)+");";
        const double gamma = 11./3.;
        //const double gamma = 2.+3./2.;
        //std::cerr<<"Kraichnan turbulence, gamma:"<<gamma<<std::endl;
        const double gamma_large_scale = -2;
        std::vector<double> norms;
        std::vector<double> Ks;
        double sumNorm = 0;
//        std::cerr<<"DEBUG BUILD: see 303 line at Defletion3D.cpp"<<std::endl;
//        #warning DEBUG BUILD: k*fLc < 1000 -> k*fLc < 100
        uniqHash += "gamma=";
        uniqHash += std::to_string(gamma);
        for(double k = 1./fLc; k*fLc < 100; k*=aMultK){ // L_max = 2п*fLc
            Ks.push_back(k);
            double normK = 0;
            normK = k*k*k*pow(k*fLc,-gamma); // Spectrum definition
            norms.push_back(normK);
            sumNorm += normK;
            fWaves.push_back(0);
        }
        //std::cerr<<"Number of waves in TurbulentMF:"<<norms.size()<<std::endl;
        uniqHash += ";Nwaves=";
        uniqHash += std::to_string(norms.size());
        for(uint i=0; i<norms.size(); i++){
//        	if(i < norms.size()-2)norms[i]=0.;
//        	else std::cerr<<"DEBUG RUN! Magnetic field was broken! Wave amp:"<<aVariance_Gauss*sqrt(norms[i]/sumNorm)<<std::endl;
            fWaves[i] = new MonochromaticMF(aRandomizer, 2.*M_PI/Ks[i]/units.Mpc, aVariance_Gauss*sqrt(norms[i]/sumNorm));
        }
        uniqHash += ";aRandomizer.Rand()=";
        uniqHash += std::to_string(aRandomizer.Rand());
        //std::cerr<<uniqHash<<std::endl;

//void test_fftw(){
//        const int N = 8;
//        fftw_complex *in, *out;
//        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//
//        fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//
//        for (int i = 0; i < N; ++i) {
//                in[i][0] = sin(2.0 * M_PI * i / N); // Вещественная часть
//                in[i][1] = 0.0;                     // Мнимая часть
//        }
//
//        fftw_execute(plan);
//
//        std::cout << "FFTW TESTING:\nАмплитуды:" << std::endl;
//        for (int i = 0; i < N/2 + 1; ++i) {
//                double amplitude = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]) / N * 2;
//                std::cout << "Частота " << i << ": " << amplitude << std::endl;
//        }
//
//        fftw_destroy_plan(plan);
//        fftw_free(in);
//        fftw_free(out);
//
//        fftw_cleanup();
//}

    }

    FFT_TurbulentMF::FFT_TurbulentMF(Randomizer &aRandomizer, double aLcor_Mpc, double aVariance_Gauss,
                    double cube_size, double aMinVarScale, int num_threads, double fSourceX, double fSourceY, double fSourceZ, FillingType mFillingT)
    {
        size_t N = (size_t)ceil(2*cube_size/(aMinVarScale*aLcor_Mpc));//2*cube_size/(aMinVarScale*aLcor_Mpc);
        if (N % 2 != 0) N++;
        ASSERT(N>2);
        ASSERT(num_threads>0);
        FillingT = mFillingT;
        uniqHash = "FFT_TurbulentMF(";
        uniqHash += std::to_string(aLcor_Mpc)+",";
        uniqHash += std::to_string(aVariance_Gauss)+",";
        uniqHash += std::to_string(cube_size)+",";
        uniqHash += std::to_string(aMinVarScale)+",";
        uniqHash += std::to_string(num_threads)+",";
        uniqHash += std::to_string(fSourceX)+",";
        uniqHash += std::to_string(fSourceY)+",";
        uniqHash += std::to_string(fSourceZ)+",";
        uniqHash += std::to_string(FillingT)+");N=";
        uniqHash += std::to_string(N);
        std::cerr<<"N = "<<N<<std::endl;

        x_min = (fSourceX-cube_size)*units.Mpc;
        y_min = (fSourceY-cube_size)*units.Mpc;
        z_min = (fSourceZ-cube_size)*units.Mpc;
        x_max = (fSourceX+cube_size)*units.Mpc;
        y_max = (fSourceY+cube_size)*units.Mpc;
        z_max = (fSourceZ+cube_size)*units.Mpc;
        Nx = N;
        Ny = N;
        Nz = N;

        int init_result = fftw_init_threads();

        if(init_result == 0){
        	std::cerr<<"FFTW threading failed!"<<std::endl;
        	exit(1);
        }
		omp_set_num_threads(num_threads);
		fftw_plan_with_nthreads(num_threads);

        if(!fftw_import_wisdom_from_filename("fftw_wisdom.fftw")){
        	std::cerr<<"Failed to load FFTW wisdom. Ignoring."<<std::endl;
        };

		int time_limit = (std::log(N)+5)*10;

		fftw_set_timelimit(time_limit);

        fftw_complex *Bx_k, *By_k, *Bz_k;
        fftw_plan plan_Bx, plan_By, plan_Bz;

        size_t complex_size = N * N * (floor(N/2) + 1);
        size_t inplace_size = N * N * (N + 2);
//        Bx_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * complex_size);
//        By_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * complex_size);
//        Bz_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * complex_size);
		std::cerr<<"Allocating memory..."<<std::endl;
        Bx_k = (fftw_complex*) fftw_malloc(sizeof(double) * inplace_size*3);
        if(!Bx_k){
        	std::cerr<<"Cannot allocate memory! Exit."<<std::endl;
        	exit(2);
        }
        By_k = Bx_k+complex_size;//(fftw_complex*) fftw_malloc(sizeof(double) * inplace_size);
        Bz_k = By_k+complex_size;//(fftw_complex*) fftw_malloc(sizeof(double) * inplace_size);

        size_t real_size = N * N * N;

        B_X = (double*) Bx_k;//(double*) fftw_malloc(sizeof(double) * real_size *3);
        B_Y = (double*) By_k;//B_X+real_size;//(double*) By_k;//fftw_malloc(sizeof(double) * real_size);
        B_Z = (double*) Bz_k;//B_Y+real_size;//(double*) Bz_k;//fftw_malloc(sizeof(double) * real_size);

        std::cerr<<"sizeof(fftw_complex) * complex_size = "<<sizeof(fftw_complex) * complex_size<<std::endl;
        std::cerr<<"sizeof(double) * real_size = "<<sizeof(double) * real_size<<std::endl;
        std::cerr<<"sizeof(double) * inplace_size = "<<sizeof(double) * inplace_size<<std::endl;
        std::cerr<<"Planning FFTW (no more "<<time_limit<<" seconds)..."<<std::endl;

        plan_Bx = fftw_plan_dft_c2r_3d(N, N, N, Bx_k, B_X, FFTW_ESTIMATE);
        plan_By = fftw_plan_dft_c2r_3d(N, N, N, By_k, B_Y, FFTW_ESTIMATE);
        plan_Bz = fftw_plan_dft_c2r_3d(N, N, N, Bz_k, B_Z, FFTW_ESTIMATE);

//		int dims[3] = {(int)N, (int)N, (int)N};
//		int inembed[3] = {(int)N, (int)N, (int)(N/2 + 1)}; // Последняя размерность комплексного массива "обрезана"
//		int onembed[3] = {(int)N, (int)N, (int)N};
//        fftw_plan plan = fftw_plan_many_dft_c2r(
//			3, dims, 3,
//			Bx_k, inembed, 1, complex_size,
//			B_X, onembed, 1, inplace_size,
//			FFTW_ESTIMATE
//		);

//		fftw_iodim64 dims[3];
//		dims[0].n = N;
//		dims[0].is = N * (N/2 + 1);   // шаг по первому измерению для комплексного массива
//		dims[0].os = N * N;           // шаг по первому измерению для реального массива
//
//		dims[1].n = N;
//		dims[1].is = (N/2 + 1);       // шаг по второму измерению для комплексного массива
//		dims[1].os = N;               // шаг по второму измерению для реального массива
//
//		dims[2].n = N;
//		dims[2].is = 1;               // шаг по последнему измерению для комплексного массива
//		dims[2].os = 1;               // шаг по последнему измерению для реального массива
//		fftw_iodim64 howmany_dims[1];
//		howmany_dims[0].n = 3;
//		howmany_dims[0].is = complex_size;   // расстояние между комплексными массивами
//		howmany_dims[0].os = inplace_size;   // расстояние между реальными массивами
//		fftw_plan plan = fftw_plan_guru64_dft_c2r(
//			3, dims, 1, howmany_dims,
//			Bx_k, B_X,
//			FFTW_MEASURE
//		);
//		if (!plan) {
//		    std::cerr << "FFTW plan creation failed!" << std::endl;
//		    exit(3);
//		}
        if(!fftw_export_wisdom_to_filename("fftw_wisdom.fftw")){
        	std::cerr<<"N = "<<N<<"; num_threads = "<<num_threads<<std::endl;
        	std::cerr<<"Failed to export FFTW wisdom. Ignoring."<<std::endl;
        };

        std::memset(Bx_k, 0, sizeof(double) * inplace_size *3); // set zeros
        //std::memset(By_k, 0, sizeof(double) * inplace_size);
        //std::memset(Bz_k, 0, sizeof(double) * inplace_size);

        double k_min = 2.0 * M_PI / (2*cube_size);
        double k_max = M_PI / aMinVarScale;
        std::cerr<<"k_min = "<<k_min<<std::endl;
        std::cerr<<"k_max = "<<k_max<<std::endl;

        const double gamma = 11./3.;
        double normalizer = 1/pow(2*M_PI,3);

        unsigned long int* aRandomizerS = (unsigned long int*)malloc(sizeof(unsigned long int)*N);

        for(int i=0;i<N;i++)aRandomizerS[i] = aRandomizer.CreateIndependent();

        uniqHash += ";k_min=";
        uniqHash += std::to_string(k_min);
        uniqHash += ";k_max=";
        uniqHash += std::to_string(k_max);
        uniqHash += ";aRandomizerS[0]=";
        uniqHash += std::to_string(aRandomizerS[0]);
        uniqHash += ";aRandomizerS[N-1]=";
        uniqHash += std::to_string(aRandomizerS[N-1]);
        std::cerr<<"Spectrum building..."<<std::endl;

		#pragma omp parallel for
        for (size_t i = 0; i < N; ++i){
        	Randomizer cur_randomizer(aRandomizerS[i]);
   	        if(omp_get_thread_num()==0 && i==0){
   	        	std::cerr<<"omp_get_num_threads() -> "<<omp_get_num_threads()<<std::endl;
   	        	std::cerr<<"omp_get_max_threads() -> "<<omp_get_max_threads()<<std::endl;
   	        }
			long double kx = (i <= N/2) ? k_min*i : -k_min*(N-i);
        	for (size_t j = 0; j < N; ++j) {
				long double ky = (j <= N/2) ? k_min*j : -k_min*(N-j);
        		for (size_t k = 0; k < (floor(N/2) + 1); ++k) {
        			long double kz = k_min*k;
        			size_t ind3d = k + (N/2+1) * (j + N * i);//ijk2ind( i , j , k );
        			long double k_mag = std::sqrt(kx*kx + ky*ky + kz*kz); // Physical k
//        			if (k_mag < k_min || k_mag > k_max) continue;
        			if (k_mag < k_min) continue;
//        			if (k_mag < 1./aLcor_Mpc) continue;
        			if (k_mag < 1./aLcor_Mpc || k_mag > 100./aLcor_Mpc) continue;
        			// Spectrum is here!
        			long double amplitude = k_mag*k_mag*pow(k_mag/aLcor_Mpc,-gamma);// * aVariance_Gauss * normalizer;
        			// divB=0 => kx*Bx(k)+ky*By(k)+kz*Bz(k)=0
        			long double nx = kx / k_mag;
        			long double ny = ky / k_mag;
        			long double nz = kz / k_mag;
      				long double tmpB[3];
        			long double tmp_norm=0;
					long double Bv[3] = {0,0,0};
					do{
						for(size_t Bvind=0;Bvind<3;Bvind++){
							Bv[Bvind] = 1.-cur_randomizer.Rand()*2.;
						}
						tmp_norm = sqrt(Bv[0]*Bv[0]+Bv[1]*Bv[1]+Bv[2]*Bv[2]);
						Bv[0] /= tmp_norm;
						Bv[1] /= tmp_norm;
						Bv[2] /= tmp_norm;
						tmp_norm = nx*Bv[0]+ny*Bv[1]+nz*Bv[2];
					}while(1-abs(tmp_norm)<1e-3); // Too parallel
					Bv[0] -= tmp_norm*nx;
					Bv[1] -= tmp_norm*ny;
					Bv[2] -= tmp_norm*nz;
					tmp_norm = sqrt(Bv[0]*Bv[0]+Bv[1]*Bv[1]+Bv[2]*Bv[2]);
					Bv[0] /= tmp_norm;
					Bv[1] /= tmp_norm;
					Bv[2] /= tmp_norm;
					if( std::isnan(Bv[0]) ){
						std::cerr<<"NaN."<<std::endl;
						std::cerr<<" tmp_norm = "<<tmp_norm<<std::endl;
						exit(1);
					}
					if( abs(nx*Bv[0]+ny*Bv[1]+nz*Bv[2])>1e-16 ){
        				std::cerr<<"..."<<std::endl;
        				std::cerr<<" tmp_norm = "<<tmp_norm<<std::endl;
        				std::cerr<<" nx*Bv = "<<nx*Bv[0]+ny*Bv[1]+nz*Bv[2]<<std::endl;
        				std::cerr<<"n  = ["<<nx<<","<<ny<<","<<nz<<"]"<<std::endl;
        				std::cerr<<"Bv = ["<<Bv[0]<<","<<Bv[1]<<","<<Bv[2]<<"]"<<std::endl;
        			}
        			long double phase2 = cur_randomizer.Rand()*M_PI*2.;
					for(size_t Bvind=0;Bvind<3;Bvind++){
						tmpB[Bvind] = amplitude*Bv[Bvind];
					}
        			Bx_k[ ind3d ][0] = tmpB[0]*cos(phase2); // Real part
        			By_k[ ind3d ][0] = tmpB[1]*cos(phase2);
        			Bz_k[ ind3d ][0] = tmpB[2]*cos(phase2);
        			Bx_k[ ind3d ][1] = tmpB[0]*sin(phase2); // Imaginary part
        			By_k[ ind3d ][1] = tmpB[1]*sin(phase2);
        			Bz_k[ ind3d ][1] = tmpB[2]*sin(phase2);
        			if( std::isnan(tmpB[0]) )std::cerr<<kx<<" "<<ky<<" "<<kz<<" "<<tmpB[0]<<" "<<tmpB[1]<<" "<<tmpB[2]<<std::endl;
        		}
        	}
        }
        free(aRandomizerS);
        Bx_k[0][0] = 0.0; Bx_k[0][1] = 0.0;
        By_k[0][0] = 0.0; By_k[0][1] = 0.0;
        Bz_k[0][0] = 0.0; Bz_k[0][1] = 0.0;

//        std::cerr<<"Executing FFT (all)..."<<std::endl;
//        fftw_execute(plan);
//        fftw_destroy_plan(plan);

        std::cerr<<"Executing FFT (1/3)..."<<std::endl;
        fftw_execute(plan_Bx);
        fftw_destroy_plan(plan_Bx);
        //fftw_free(Bx_k);
        std::cerr<<"Executing FFT (2/3)..."<<std::endl;
        fftw_execute(plan_By);
        fftw_destroy_plan(plan_By);
        //fftw_free(By_k);
        std::cerr<<"Executing FFT (3/3)..."<<std::endl;
        fftw_execute(plan_Bz);
        fftw_destroy_plan(plan_Bz);
        //fftw_free(Bz_k);
        std::cerr<<"Executing FFT done!"<<std::endl;

        fftw_cleanup();
        fftw_cleanup_threads();

		double meanB = 0.;
		#pragma omp parallel reduction(+:meanB)
		for(size_t i=0;i<N;++i){
			for(size_t j=0;j<N;++j){
				#pragma omp for
				for(size_t k=0;k<N;++k){
					size_t ind3d_old = k + (N+2) * (j + N * i);// in-place for 3 separated plans
					size_t ind3d_new = k + N * (j + N * i);
					B_X[ind3d_new] = B_X[ind3d_old];
					B_Y[ind3d_new] = B_Y[ind3d_old];
					B_Z[ind3d_new] = B_Z[ind3d_old];
					meanB += (B_X[ind3d_new]*B_X[ind3d_new]+B_Y[ind3d_new]*B_Y[ind3d_new]+B_Z[ind3d_new]*B_Z[ind3d_new]);
				}
			}
		}
        meanB /= real_size;
        meanB = sqrt(meanB);
		std::cerr<<"meanB = "<<meanB<<".\t";
		double renorm = aVariance_Gauss/meanB;
		std::cerr<<"renorm = "<<renorm<<". Normalize..."<<std::endl;
		#pragma omp parallel for
        for(size_t i=0;i<real_size;++i){
        	B_X[i] *= renorm;
        	B_Y[i] *= renorm;
        	B_Z[i] *= renorm;
        }
    }
    FFT_TurbulentMF::~FFT_TurbulentMF(){
    	std::cerr<<"Call FFT_TurbulentMF destructor"<<std::endl;
    	fftw_free((fftw_complex*)B_X);
    	B_X = nullptr;
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
        SmartPtr<MagneticField> mf = new TurbulentMF(r, lambdaMpc, Bgauss);
        Deflection3D defl(mf);
        std::ofstream outB("TurbulentMF_B.txt");
        mf->Print(lambdaMpc,6.*lambdaMpc, outB);
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

    CellMF::CellMF(Randomizer &aRandomizer, double aLcor_Mpc, double aAmplitudeGauss):
        fLc(aLcor_Mpc*units.Mpc),
        fAmplitude(aAmplitudeGauss),
        fCurrentMaxz(0),
        fSize(10000)
    {
        double val = 0;
        for (int i = 0; i < fSize; i++)
        {
            val = 1 - aRandomizer.Rand()*2;
            if (val < 0)
                fCellDirection.push_back(-1);
            else
                fCellDirection.push_back(1);
            fCurrentMaxz += fLc;
        }
    }
    double CellMF::MinVariabilityScale(const CosmoTime &aTime) const
    {
        return 0.05*fLc;
    }
    void CellMF::GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        outValue[0] = 0;
        outValue[1] = fAmplitude*fCellDirection[int(fSize*x[2]/fCurrentMaxz)];
        outValue[2] = 0;
        // std::cout << fCurrentMaxz << '\t' << fLc << '\t' << x[2] << '\t' << int(fSize*x[2]/fCurrentMaxz) << '\t' << outValue[1] << '\n';
        // std::cin.get();
    }

//////// MF from file
    SavedMF::SavedMF(const char* fname,const FillingType mFillingType)
    {
    	origfname = fname;
    	FillingT = mFillingType;
    	std::cout << "Analyse file..." << std::endl;
        std::ifstream fl(fname);
		std::set<double> x_v;
		std::set<double> y_v;
		std::set<double> z_v;
		double tmpX,tmpY,tmpZ,tmpBx,tmpBy,tmpBz;
		std::string line;
		std::vector<std::string> lines;
		x_min = 0;
		y_min = 0;
		z_min = 0;
		x_max = 0;
		y_max = 0;
		z_max = 0;
		while (std::getline(fl, line)) {
			if(line[0]=='#')continue; // comment line
			lines.push_back(std::string(line));
			std::istringstream iss(line);
			iss >> tmpX >> tmpY >> tmpZ >> tmpBx >> tmpBy >> tmpBz;
			x_v.insert(tmpX);
			y_v.insert(tmpY);
			z_v.insert(tmpZ);
			if(tmpX > x_max)x_max = tmpX;
			else if(tmpX < x_min)x_min = tmpX;
			if(tmpY > y_max)y_max = tmpY;
			else if(tmpY < y_min)y_min = tmpY;
			if(tmpZ > z_max)z_max = tmpZ;
			else if(tmpZ < z_min)z_min = tmpZ;
        }
        fl.close();
        //fl.clear();
		//fl.seekg(0);
        Nx = x_v.size();
        Ny = y_v.size();
        Nz = z_v.size();
        std::cout << " Nx = " << Nx;
        std::cout << "; Ny = " << Ny;
        std::cout << "; Nz = " << Nz << std::endl;
        std::cout << " x_min = " << x_min;
        std::cout << "; y_min = " << y_min;
        std::cout << "; z_min = " << z_min << std::endl;
        std::cout << " x_max = " << x_max;
        std::cout << "; y_max = " << y_max;
        std::cout << "; z_max = " << z_max << std::endl;
        x_min*=units.Mpc;
        y_min*=units.Mpc;
        z_min*=units.Mpc;
        x_max*=units.Mpc;
        y_max*=units.Mpc;
        z_max*=units.Mpc;
        B_X = new double[Nx*Ny*Nz];
        B_Y = new double[Nx*Ny*Nz];
        B_Z = new double[Nx*Ny*Nz];
        std::cout << "Reading field..." << std::endl;
        for(size_t ii=0;ii<lines.size();ii++){
        	line = lines[ii];
		//while (std::getline(fl, line)) {
			if(line[0]=='#')continue; // comment line
			std::istringstream iss(line);
			iss >> tmpX >> tmpY >> tmpZ >> tmpBx >> tmpBy >> tmpBz;
			size_t ind = coord2index(tmpX*units.Mpc, tmpY*units.Mpc, tmpZ*units.Mpc);
			B_X[ind] = tmpBx;
			B_Y[ind] = tmpBy;
			B_Z[ind] = tmpBz;
        }
        //fl.close();
        std::cout << "Done!" << std::endl;
    }
    size_t SavedMF::ijk2ind(size_t i, size_t j, size_t k) const
    {
    	return k + Nz * (j + Ny * i);
    }
    size_t SavedMF::coord2index(double x, double y, double z) const
    {
    	int i,j,k;
    	if(x<=x_min)i = 0;
    	else if(x>=x_max)i = Nx-1;
    	else i = round( (double(Nx-1))*(x-x_min)/(x_max-x_min) );

    	if(y<=y_min)j = 0;
    	else if(y>=y_max)j = Ny-1;
    	else j = round( (double(Ny-1))*(y-y_min)/(y_max-y_min) );

    	if(z<=z_min)k = 0;
    	else if(z>=z_max)k = Nz-1;
    	else k = round( (double(Nz-1))*(z-z_min)/(z_max-z_min) );
    	return ijk2ind(i,j,k);
//    	return k + Nz * (j + Ny * i);
//        return i+Nx*(j+Ny*k); // (Nx-1)+Nx*(Ny-1)+Nx*Ny*(Nz-1) = Nx*Ny*Nz-1
    }
    size_t SavedMF::coord2indexZero(double coord, double min_val, double max_val, size_t N) const
    {
    	if(coord<=min_val)return 0;
    	else if(coord>=max_val)return N-1;
    	else return round( (double(N-1))*(coord-min_val)/(max_val-min_val) );
    }
    void SavedMF::GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        ASSERT(outValue.size()==3);
        GetValueGauss_LINEAR(x, aTime, outValue);
        //GetValueGauss_NEAREST(x, aTime, outValue);
    }
    void SavedMF::GetValueGauss_NEAREST(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        size_t ind = coord2index(x[0], x[1], x[2]);
        outValue[0] = B_X[ind];
        outValue[1] = B_Y[ind];
        outValue[2] = B_Z[ind];
    }
    double SavedMF::linApprox(const double* f,
    					size_t i, size_t j, size_t k,
    					double x, double y, double z) const
    {
		if (i >= Nx-1) i = Nx-2;
		if (j >= Ny-1) j = Ny-2;
		if (k >= Nz-1) k = Nz-2;
		double x0 = x_min+double(i)*(x_max-x_min)/double(Nx-1);
		double x1 = x_min+double(i+1)*(x_max-x_min)/double(Nx-1);
		double y0 = y_min+double(j)*(y_max-y_min)/double(Ny-1);
		double y1 = y_min+double(j+1)*(y_max-y_min)/double(Ny-1);
		double z0 = z_min+double(k)*(z_max-z_min)/double(Nz-1);
		double z1 = z_min+double(k+1)*(z_max-z_min)/double(Nz-1);
//    	double C000 = f[ i  +Nx*(j  +Ny*  k  ) ];
//    	double C100 = f[ i+1+Nx*(j  +Ny*  k  ) ];
//    	double C010 = f[ i  +Nx*(j+1+Ny*  k  ) ];
//    	double C001 = f[ i  +Nx*(j  +Ny*(k+1)) ];
//    	double C110 = f[ i+1+Nx*(j+1+Ny*  k  ) ];
//    	double C011 = f[ i  +Nx*(j+1+Ny*(k+1)) ];
//    	double C101 = f[ i+1+Nx*(j  +Ny*(k+1)) ];
//    	double C111 = f[ i+1+Nx*(j+1+Ny*(k+1)) ];
//    	double C000 = f[ k  +Nz*(j  +Ny*  i  ) ];
//    	double C100 = f[ k  +Nz*(j  +Ny*(i+1)) ];
//    	double C010 = f[ k  +Nz*(j+1+Ny*  i  ) ];
//    	double C001 = f[ k+1+Nz*(j  +Ny*  i  ) ];
//    	double C110 = f[ k  +Nz*(j+1+Ny*(i+1)) ];
//    	double C011 = f[ k+1+Nz*(j+1+Ny*  i  ) ];
//    	double C101 = f[ k+1+Nz*(j  +Ny*(i+1)) ];
//    	double C111 = f[ k+1+Nz*(j+1+Ny*(i+1)) ];
    	double C000 = f[ ijk2ind( i , j , k ) ];
    	double C100 = f[ ijk2ind(i+1, j , k ) ];
    	double C010 = f[ ijk2ind( i ,j+1, k ) ];
    	double C001 = f[ ijk2ind( i , j ,k+1) ];
    	double C110 = f[ ijk2ind(i+1,j+1, k ) ];
    	double C011 = f[ ijk2ind( i ,j+1,k+1) ];
    	double C101 = f[ ijk2ind(i+1, j ,k+1) ];
    	double C111 = f[ ijk2ind(i+1,j+1,k+1) ];
		double denom = (x0-x1)*(y0-y1)*(z0-z1);
		double a0 = -C000*x1*y1*z1 + C001*x1*y1*z0 + C010*x1*y0*z1 - C011*x1*y0*z0
					+C100*x0*y1*z1 - C101*x0*y1*z0 - C110*x0*y0*z1 + C111*x0*y0*z0;
		double a1 = C000*y1*z1 - C001*y1*z0 - C010*y0*z1 + C011*y0*z0
					-C100*y1*z1 + C101*y1*z0 + C110*y0*z1 - C111*y0*z0;
		double a2 = C000*x1*z1 - C001*x1*z0 - C010*x1*z1 + C011*x1*z0
					-C100*x0*z1 + C101*x0*z0 + C110*x0*z1 - C111*x0*z0;
		double a3 = C000*x1*y1 - C001*x1*y1 - C010*x1*y0 + C011*x1*y0
					-C100*x0*y1 + C101*x0*y1 + C110*x0*y0 - C111*x0*y0;
		double a4 = -C000*z1 + C001*z0 + C010*z1 - C011*z0
					+C100*z1 - C101*z0 - C110*z1 + C111*z0;
		double a5 = -C000*y1 + C001*y1 + C010*y0 - C011*y0
					+C100*y1 - C101*y1 - C110*y0 + C111*y0;
		double a6 = -C000*x1 + C001*x1 + C010*x1 - C011*x1
					+C100*x0 - C101*x0 - C110*x0 + C111*x0;
		double a7 = C000 - C001 - C010 + C011
					-C100 + C101 + C110 - C111;
		return (a0 + a1*x + a2*y + a3*z
				+ a4*x*y + a5*x*z + a6*y*z + a7*x*y*z)/denom;
    }
    void SavedMF::GetValueGauss_OUTOFBOUND(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        //std::cerr<<"Warning: out of bound fallback: ("<<x[0]/units.kpc<<","<<x[1]/units.kpc<<","<<x[2]/units.kpc<<") kpc;"<<std::endl;
        switch (FillingT) {
        	case FILLING_ZERO:
        		outValue[0] = 0;
				outValue[1] = 0;
				outValue[2] = 0;
        		break;
        	case FILLING_NEAREST:
        		GetValueGauss_NEAREST(x, aTime, outValue);
        		break;
        	case FILLING_REPEAT:
        		{
        		long double* newX = new long double[3];
        		double dx = x_max-x_min;
        		double dy = y_max-y_min;
        		double dz = z_max-z_min;
//        		newX[0] = x[0];
//        		newX[1] = x[1];
//        		newX[2] = x[2];
				long double diffX = x[0] - x_min;
				long double remX = fmodl(diffX, dx);
				if (remX < 0) remX += dx;
				newX[0] = x_min + remX;

				// Приведение y
				long double diffY = x[1] - y_min;
				long double remY = fmodl(diffY, dy);
				if (remY < 0) remY += dy;
				newX[1] = y_min + remY;

				// Приведение z
				long double diffZ = x[2] - z_min;
				long double remZ = fmodl(diffZ, dz);
				if (remZ < 0) remZ += dz;
				newX[2] = z_min + remZ;
//        		newX[0] = x[0]-dx*std::round(x[0]/dx);
//        		newX[1] = x[1]-dy*std::round(x[1]/dy);
//        		newX[2] = x[2]-dz*std::round(x[2]/dz);
//				while(newX[0]>=x_max)newX[0]-=dx;// TODO: Optimize it
//				while(newX[0]<=x_min)newX[0]+=dx;
//				while(newX[1]>=y_max)newX[1]-=dy;// TODO: Optimize it
//				while(newX[1]<=y_min)newX[1]+=dy;
//				while(newX[2]>=z_max)newX[2]-=dz;// TODO: Optimize it
//				while(newX[2]<=z_min)newX[2]+=dz;
        		// На стыке кубов будет лажа, но, видимо, это неизбежность
        		if(newX[0]>=x_max)newX[0]-=dx/Nx;
        		else if(newX[0]<=x_min)newX[0]+=dx/Nx;
        		if(newX[1]>=y_max)newX[1]-=dy/Ny;
        		else if(newX[1]<=y_min)newX[1]+=dy/Ny;
        		if(newX[2]>=z_max)newX[2]-=dz/Nz;
        		else if(newX[2]<=z_min)newX[2]+=dz/Nz;
        		GetValueGauss_LINEAR(newX, aTime, outValue);
        		delete [] newX;
        		}
        		break;
        	case FILLING_MIRROR:
        		{
        		long double* newX = new long double[3];
        		double dx = x_max-x_min;
        		double dy = y_max-y_min;
        		double dz = z_max-z_min;
				long double diffX = x[0] - x_min;
				long double remX = fmodl(diffX, 2.*dx);
				if (remX < 0) remX += 2.*dx;
				if(remX<dx)newX[0] = x_min + remX;
				else newX[0] = x_min + (2.*dx-remX);

				// Приведение y
				long double diffY = x[1] - y_min;
				long double remY = fmodl(diffY, 2.*dy);
				if (remY < 0) remY += 2.*dy;
				if(remY<dy)newX[1] = y_min + remY;
				else newX[1] = y_min + (2.*dy-remY);

				// Приведение z
				long double diffZ = x[2] - z_min;
				long double remZ = fmodl(diffZ, 2.*dz);
				if (remZ < 0) remZ += 2.*dz;
				if(remZ<dz)newX[2] = z_min + remZ;
				else newX[2] = z_min + (2.*dz-remZ);
        		// На стыке кубов будет лажа, но, видимо, это неизбежность
        		if(newX[0]>=x_max)newX[0]-=dx/Nx;
        		else if(newX[0]<=x_min)newX[0]+=dx/Nx;
        		if(newX[1]>=y_max)newX[1]-=dy/Ny;
        		else if(newX[1]<=y_min)newX[1]+=dy/Ny;
        		if(newX[2]>=z_max)newX[2]-=dz/Nz;
        		else if(newX[2]<=z_min)newX[2]+=dz/Nz;
        		GetValueGauss_LINEAR(newX, aTime, outValue);
        		delete [] newX;
        		}
        		break;
        }
		return;
	}
    void SavedMF::GetValueGauss_LINEAR(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
    	//I need all i,j,k coord
    	size_t i,j,k;
    	if(x[0]<=x_min){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;} // Use nearest out of the bound
    	else if(x[0]>=x_max){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;}
    	else i = size_t( (double(Nx-1))*(x[0]-x_min)/(x_max-x_min) );

    	if(x[1]<=y_min){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;}
    	else if(x[1]>=y_max){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;}
    	else j = size_t( (double(Ny-1))*(x[1]-y_min)/(y_max-y_min) );

    	if(x[2]<=z_min){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;}
    	else if(x[2]>=z_max){GetValueGauss_OUTOFBOUND(x, aTime, outValue);return;}
    	else k = size_t( (double(Nz-1))*(x[2]-z_min)/(z_max-z_min) );

    	// Use nearest near the bound
    	//if( i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1 ){GetValueGauss_NEAREST(x, aTime, outValue);return;}
        outValue[0] = linApprox(B_X, i,j,k, x[0],x[1],x[2]);
        outValue[1] = linApprox(B_Y, i,j,k, x[0],x[1],x[2]);
        outValue[2] = linApprox(B_Z, i,j,k, x[0],x[1],x[2]);
    }
    double SavedMF::MinVariabilityScale(const CosmoTime &aTime) const
    {
    	double a=(x_max-x_min)/Nx;
    	double b=(y_max-y_min)/Ny;
    	double c=(z_max-z_min)/Nz;
    	if (a<b && a<c) return a;
    	if (b<a && b<c) return b;
    	return c;
    }
    double SavedMF::MinVariabilityScale(const long double *x, const CosmoTime &aTime,const long double eps) const
    {
    	double tmp = MinVariabilityScale(aTime);
    	if(FillingT == FILLING_ZERO){
    		if(x[0]<x_min-tmp-eps || x[0]>x_max+tmp+eps)return 0;
    		if(x[1]<y_min-tmp-eps || x[1]>y_max+tmp+eps)return 0;
    		if(x[2]<z_min-tmp-eps || x[2]>z_max+tmp+eps)return 0;
    	}
    	return tmp;
    }
    SavedMF::~SavedMF()
    {
    	std::cerr<<"Call SavedMF destructor"<<std::endl;
    	if(B_X==nullptr)return; // It was free before
        delete [] B_X;
        delete [] B_Y;
        delete [] B_Z;
    }

    void ComboMF::GetValueGauss(const long double *x, const CosmoTime &aTime,
                                   std::vector<double> &outValue) const
    {
         for(size_t it=0;it<3;it++)outValue[it]=0;
         std::vector<double> tmpValue(3, 0.0);
         for(size_t i = 0; i < mfv.size(); ++i){
         	mfv[i]->GetValueGauss(x,aTime,tmpValue);
         	for(size_t it=0;it<3;it++)outValue[it]+=tmpValue[it];
         }
    }
    double ComboMF::MinVariabilityScale(const CosmoTime &aTime) const
    {
         double mMinVariabilityScale=mfv[0]->MinVariabilityScale(aTime);
         double tmp;
         for(size_t i = 1; i < mfv.size(); ++i){
         	tmp = mfv[i]->MinVariabilityScale(aTime);
         	if(tmp>0 && tmp<mMinVariabilityScale)mMinVariabilityScale=tmp;
         }
         return mMinVariabilityScale;
    }
    double ComboMF::MinVariabilityScale(const long double *x,const CosmoTime &aTime,const long double eps) const
    {
         double mMinVariabilityScale=mfv[0]->MinVariabilityScale(x,aTime,eps);
         double tmp;
         for(size_t i = 1; i < mfv.size(); ++i){
         	tmp = mfv[i]->MinVariabilityScale(x,aTime,eps);
         	if(tmp>0 && tmp<mMinVariabilityScale)mMinVariabilityScale=tmp;
         }
         return mMinVariabilityScale;
    }
    void ComboMF::PrintUniqHash(std::ostream& aOutput) const
    {
    	aOutput<<"ComboMF{"<<std::endl;
        for(size_t i = 0; i < mfv.size(); ++i){
        	aOutput<<"\t";
        	mfv[i]->PrintUniqHash(aOutput);
        }
        aOutput<<"}"<<std::endl;
    }


//////// Galactic MF
    void GalacticMF::GetValueGauss(const long double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        outValue[0] = 0;
        outValue[1] = 0;
        outValue[2] = 0;
        // std::cout << fCurrentMaxz << '\t' << fLc << '\t' << x[2] << '\t' << int(fSize*x[2]/fCurrentMaxz) << '\t' << outValue[1] << '\n';
        // std::cin.get();
    }
    double GalacticMF::MinVariabilityScale(const CosmoTime &aTime) const
    {
        return 0.05*units.kpc; // TODO
    }


}
