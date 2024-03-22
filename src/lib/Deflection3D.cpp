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
        const coord_type* N = aParticle.Pdir;
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
        coord_type newPdir[3];
        int totSteps=0;//debug info

        for(u_int stepB = 0; stepB<nSteps; stepB++){
            pB->GetValueGauss(aParticle.X, aParticle.Time, curB);
            GetRotationRate(aParticle, curB, deflRate);
            cosmo_time absRate = sqrt(deflRate[0]*deflRate[0]+deflRate[1]*deflRate[1]+deflRate[2]*deflRate[2]);
            cosmo_time maxStep = 0.05*M_PI/absRate/fAccuracy;//max step corresponds to 9/fAccuracy degrees turn
            u_int nSubSteps = (u_int)(dx/maxStep+1.);
            cosmo_time dXsubStep = dx/nSubSteps;
            for(u_int substep=0; substep<nSubSteps; substep++){
                if(substep)
                    GetRotationRate(aParticle, curB, deflRate);
                cosmo_time absPdir = 0.;
                for(int i=0; i<3; i++){
                    coord_type val = aParticle.Pdir[i] + deflRate[i]*dXsubStep;
                    newPdir[i] = val;
                    absPdir += val*val;
                }
                absPdir = sqrt(absPdir);//norm of new vector
                for(int i=0; i<3; i++){
                    coord_type val = newPdir[i]/absPdir;//make sure rotated vector has unit norm
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

    void MonochromaticMF::GetValueGauss(const coord_type *x, const CosmoTime &aTime,
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

    int MagneticField::Print(double lambdaMpc, std::ostream& aOutput, coord_type step) {
        std::vector<double> B(3);
        CosmoTime t(0);
        if(step<0)step=lambdaMpc*units.Mpc/20.; // Default step -- 1/20
        else step*=units.Mpc;
        coord_type L=lambdaMpc*units.Mpc*3;
        aOutput << "#x,Mpc" << "\t" << "y,Mpc" << "\t" << "z,Mpc" << "\t" << "Bx" << "\t" << "By" << "\t" << "Bz" << "\n";
        for(coord_type x=-L; x<=L; x+=step) {
        	std::cerr<<"Saved " << double(int( (L+x)/(2*L) * 1000))/10. << "%...\r";
        	std::cerr.flush();
            for (coord_type y = -L; y <= L; y+=step) {
                for (coord_type z = -L; z <= L; z += step) {
                    const coord_type r[] = {x, y, z};
                    GetValueGauss(r, t, B);
                    aOutput << x / units.Mpc << "\t" << y / units.Mpc << "\t" << z / units.Mpc << "\t"
                              << B[0] << "\t" << B[1] << "\t" << B[2] << "\n";
                }
                aOutput << "\n";
            }
            aOutput << "\n\n";
        }
        aOutput << "# use gnuplot command \"set pm3d at b; splot 'magneticField.out' i 0 u 2:3:col\"  where col=4,5 or 6 to view the data\n";
        std::cerr<<"Saved 100%   "<<std::endl;
        return 0;
    }

    double MagneticField::CorrelationLength(std::vector<coord_type> x, std::vector<coord_type> aDirection, const CosmoTime &aTime,
                                            double aMaxValMpc) {
        double step = 0.1*MinVariabilityScale(aTime);
        int max_steps = (int)(aMaxValMpc*units.Mpc/step + 0.5);
        if (max_steps == 0)
            return step/units.Mpc;
        auto norm = sqrt(aDirection[0]*aDirection[0] + aDirection[1]*aDirection[1] + aDirection[2]*aDirection[2]);
        for(int dim_i=0; dim_i<3; dim_i++)
            aDirection[dim_i] *= (step/norm);

        auto iniB = std::vector<double>(3);
        auto curB = std::vector<double>(3);
        GetValueGauss(x.data(), aTime, iniB);
        for(int step_i=0; step_i<max_steps; step_i++){
            for(int dim_i=0; dim_i<3; dim_i++)
                x[dim_i] += aDirection[dim_i];
            GetValueGauss(x.data(), aTime, curB);
            double sc_product = 0.;
            for(int dim_i=0; dim_i<3; dim_i++)
                sc_product += curB[dim_i]*iniB[dim_i];
            if (sc_product < 0)
                return step*(step_i+1)/units.Mpc;
        }
        return aMaxValMpc;
    }


    int MonochromaticMF::UnitTest() {
        std::vector<double> B(3);
        CosmoTime t(0);
        Randomizer r;
        double lambdaMpc = 1.;
        double Bgauss = 10.;
        MonochromaticMF mf(r,lambdaMpc,Bgauss);
        mf.Print(lambdaMpc);
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

    void TurbulentMF::GetValueGauss(const coord_type *x, const CosmoTime &aTime, std::vector<double> &outValue) const {
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

    double TurbulentMF::MeanCorLength(Randomizer &aRandomizer, double aLmin_Mpc, double aLmax_Mpc, double aVariance_Gauss,
                                      int aNmodes, const CosmoTime &aTime, int aNsamples, double aMaxValMpc){
        double l_cor_sum = 0.;
        auto x = std::vector<coord_type>({0.,0.,0.});
        auto dx = std::vector<coord_type>({0.,0.,1.});
        for(int i=0; i<aNsamples; i++){
            TurbulentMF B(aRandomizer, aLmin_Mpc, aLmax_Mpc, aVariance_Gauss, aNmodes);
            l_cor_sum += B.CorrelationLength(x, dx, aTime, aMaxValMpc);
        }
        return l_cor_sum/aNsamples;
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

//////// MF from file
    SavedMF::SavedMF(const char* fname)
    {
    	std::cerr << "Analyse file..." << std::endl;
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
        std::cerr << " Nx = " << Nx;
        std::cerr << "; Ny = " << Ny;
        std::cerr << "; Nz = " << Nz << std::endl;
        std::cerr << " x_min = " << x_min;
        std::cerr << "; y_min = " << y_min;
        std::cerr << "; z_min = " << z_min << std::endl;
        std::cerr << " x_max = " << x_max;
        std::cerr << "; y_max = " << y_max;
        std::cerr << "; z_max = " << z_max << std::endl;
        x_min*=units.Mpc;
        y_min*=units.Mpc;
        z_min*=units.Mpc;
        x_max*=units.Mpc;
        y_max*=units.Mpc;
        z_max*=units.Mpc;
        B_X = new double[Nx*Ny*Nz];
        B_Y = new double[Nx*Ny*Nz];
        B_Z = new double[Nx*Ny*Nz];
        std::cerr << "Reading field..." << std::endl;
        for(int ii=0;ii<lines.size();ii++){
        	line = lines[ii];
		//while (std::getline(fl, line)) {
			if(line[0]=='#')continue; // comment line
			std::istringstream iss(line);
			iss >> tmpX >> tmpY >> tmpZ >> tmpBx >> tmpBy >> tmpBz;
			int ind = coord2index(tmpX*units.Mpc, tmpY*units.Mpc, tmpZ*units.Mpc);
			B_X[ind] = tmpBx;
			B_Y[ind] = tmpBy;
			B_Z[ind] = tmpBz;
        }
        //fl.close();
        std::cerr << "Done!" << std::endl;
    }
    int SavedMF::coord2index(coord_type x, coord_type y, coord_type z) const
    {
    	int i,j,k;
    	if(x<=x_min)i = 0;
    	else if(x>=x_max)i = Nx-1;
    	else i = round( (coord_type(Nx-1))*(x-x_min)/(x_max-x_min) );

    	if(y<=y_min)j = 0;
    	else if(y>=y_max)j = Ny-1;
    	else j = round( (coord_type(Ny-1))*(y-y_min)/(y_max-y_min) );

    	if(z<=z_min)k = 0;
    	else if(z>=z_max)k = Nz-1;
    	else k = round( (coord_type(Nz-1))*(z-z_min)/(z_max-z_min) );
        return i+Nx*(j+Ny*k); // (Nx-1)+Nx*(Ny-1)+Nx*Ny*(Nz-1) = Nx*Ny*Nz-1
    }
    int SavedMF::coord2indexZero(double coord, double min_val, double max_val, int N) const
    {
    	if(coord<=min_val)return 0;
    	else if(coord>=max_val)return N-1;
    	else return round( (double(N-1))*(coord-min_val)/(max_val-min_val) );
    }
    void SavedMF::GetValueGauss(const coord_type *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        ASSERT(outValue.size()==3);
        GetValueGauss_LINEAR(x, aTime, outValue);
        //GetValueGauss_NEAREST(x, aTime, outValue);
    }
    void SavedMF::GetValueGauss_NEAREST(const coord_type *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        int ind = coord2index(x[0], x[1], x[2]);
        outValue[0] = B_X[ind];
        outValue[1] = B_Y[ind];
        outValue[2] = B_Z[ind];
    }
    double SavedMF::linApprox(const double* f,
    					int i, int j, int k,
    					coord_type x, coord_type y, coord_type z) const
    {
		coord_type x0 = x_min+coord_type(i)*(x_max-x_min)/coord_type(Nx-1);
		coord_type x1 = x_min+coord_type(i+1)*(x_max-x_min)/coord_type(Nx-1);
		coord_type y0 = y_min+coord_type(j)*(y_max-y_min)/coord_type(Ny-1);
		coord_type y1 = y_min+coord_type(j+1)*(y_max-y_min)/coord_type(Ny-1);
		coord_type z0 = z_min+coord_type(k)*(z_max-z_min)/coord_type(Nz-1);
		coord_type z1 = z_min+coord_type(k+1)*(z_max-z_min)/coord_type(Nz-1);
    	coord_type C000 = f[ i  +Nx*(j  +Ny*  k  ) ];
    	coord_type C100 = f[ i+1+Nx*(j  +Ny*  k  ) ];
    	coord_type C010 = f[ i  +Nx*(j+1+Ny*  k  ) ];
    	coord_type C001 = f[ i  +Nx*(j  +Ny*(k+1)) ];
    	coord_type C110 = f[ i+1+Nx*(j+1+Ny*  k  ) ];
    	coord_type C011 = f[ i  +Nx*(j+1+Ny*(k+1)) ];
    	coord_type C101 = f[ i+1+Nx*(j  +Ny*(k+1)) ];
    	coord_type C111 = f[ i+1+Nx*(j+1+Ny*(k+1)) ];
		coord_type denom = (x0-x1)*(y0-y1)*(z0-z1);
		coord_type a0 = -C000*x1*y1*z1 + C001*x1*y1*z0 + C010*x1*y0*z1 - C011*x1*y0*z0
					+C100*x0*y1*z1 - C101*x0*y1*z0 - C110*x0*y0*z1 + C111*x0*y0*z0;
		coord_type a1 = C000*y1*z1 - C001*y1*z0 - C010*y0*z1 + C011*y0*z0
					-C100*y1*z1 + C101*y1*z0 + C110*y0*z1 - C111*y0*z0;
		coord_type a2 = C000*x1*z1 - C001*x1*z0 - C010*x1*z1 + C011*x1*z0
					-C100*x0*z1 + C101*x0*z0 + C110*x0*z1 - C111*x0*z0;
		coord_type a3 = C000*x1*y1 - C001*x1*y1 - C010*x1*y0 + C011*x1*y0
					-C100*x0*y1 + C101*x0*y1 + C110*x0*y0 - C111*x0*y0;
		coord_type a4 = -C000*z1 + C001*z0 + C010*z1 - C011*z0
					+C100*z1 - C101*z0 - C110*z1 + C111*z0;
		coord_type a5 = -C000*y1 + C001*y1 + C010*y0 - C011*y0
					+C100*y1 - C101*y1 - C110*y0 + C111*y0;
		coord_type a6 = -C000*x1 + C001*x1 + C010*x1 - C011*x1
					+C100*x0 - C101*x0 - C110*x0 + C111*x0;
		coord_type a7 = C000 - C001 - C010 + C011
					-C100 + C101 + C110 - C111;
		return (a0 + a1*x + a2*y + a3*z
				+ a4*x*y + a5*x*z + a6*y*z + a7*x*y*z)/denom;
    }
    void SavedMF::GetValueGauss_LINEAR(const coord_type *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
    	//I need all i,j,k coord
    	int i,j,k;
    	if(x[0]<=x_min){GetValueGauss_NEAREST(x, aTime, outValue);return;} // Use nearest out of the bound
    	else if(x[0]>=x_max){GetValueGauss_NEAREST(x, aTime, outValue);return;}
    	else i = int( (coord_type(Nx-1))*(x[0]-x_min)/(x_max-x_min) );

    	if(x[1]<=y_min){GetValueGauss_NEAREST(x, aTime, outValue);return;}
    	else if(x[1]>=y_max){GetValueGauss_NEAREST(x, aTime, outValue);return;}
    	else j = int( (coord_type(Ny-1))*(x[1]-y_min)/(y_max-y_min) );

    	if(x[2]<=z_min){GetValueGauss_NEAREST(x, aTime, outValue);return;}
    	else if(x[2]>=z_max){GetValueGauss_NEAREST(x, aTime, outValue);return;}
    	else k = int( (coord_type(Nz-1))*(x[2]-z_min)/(z_max-z_min) );

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
    	if (a<b && a<c) return 0.4*a;
    	if (b<a && b<c) return 0.4*b;
    	return 0.4*c;
    }
    SavedMF::~SavedMF()
    {
        delete [] B_X;
        delete [] B_Y;
        delete [] B_Z;
    }
}
