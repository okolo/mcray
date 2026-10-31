//
// Created by ok on 06/08/2015.
//

//#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
//#include <boost/numeric/ublas/storage.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/geometry/geometry.hpp>
#include "Deflection3D.h"
#include "Test.h"
#include "hdf5.h"

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
        return 0.1*fLambda;
    }

    int MagneticField::Print(double lambdaMpc, std::ostream& aOutput) {
        std::vector<double> B(3);
        CosmoTime t(0);
        double step=lambdaMpc*units.Mpc/10.;
        double L=lambdaMpc*units.Mpc*6;
        double BxMean = 0;
        double ByMean = 0;
        double BzMean = 0;
        double BxyMean = 0;
        double BxzMean = 0;
        double ByzMean = 0;
        for(double x=0.; x<=L; x+=step) {
            for (double y = 0.; y <= L; y+=step) {
                for (double z = 0.; z <= L; z += step) {
                    const double r[] = {x, y, z};
                    GetValueGauss(r, t, B);
                    // aOutput << x / units.Mpc << "\t" << y / units.Mpc << "\t" << z / units.Mpc << "\t"
                    //           << B[0] << "\t" << B[1] << "\t" << B[2] << "\n";
                    BxMean += fabs(B[0]);
                    ByMean += fabs(B[1]);
                    BzMean += fabs(B[2]);
                    BxyMean += sqrt(B[0]*B[0] + B[1]*B[1]);
                    BxzMean += sqrt(B[0]*B[0] + B[2]*B[2]);
                    ByzMean += sqrt(B[1]*B[1] + B[2]*B[2]);
                }
                // aOutput << "\n";
            }
            // aOutput << "\n\n";
        }
        BxMean /= pow(L/step + 1, 3);
        ByMean /= pow(L/step + 1, 3);
        BzMean /= pow(L/step + 1, 3);
        BxyMean /= pow(L/step + 1, 3);
        BxzMean /= pow(L/step + 1, 3);
        ByzMean /= pow(L/step + 1, 3);
        std::cout << "Mean |Bx| =  " << BxMean << '\n';
        std::cout << "Mean |By| =  " << ByMean << '\n';
        std::cout << "Mean |Bz| =  " << BzMean << '\n';
        std::cout << "Mean |Bxy| =  " << BxyMean << '\n';
        std::cout << "Mean |Bxz| =  " << BxzMean << '\n';
        std::cout << "Mean |Byz| =  " << ByzMean << '\n';
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
                outValue[i] += (wave[i]);
        }
    }

    double TurbulentMF::MinVariabilityScale(const CosmoTime &aTime) const {
        return fLc*fMinVarScale;
    }

    TurbulentMF::TurbulentMF(Randomizer &aRandomizer, double aLcor_Mpc, double aVariance_Gauss, double aMultK, double aMinVarScale):
            fLc(aLcor_Mpc*units.Mpc),
            fMinVarScale(aMinVarScale)
    {
        ASSERT(aMultK>1.);
        const double gamma = 11./3.;
        const double gamma_large_scale = -2;
        std::vector<double> norms;
        std::vector<double> Ks;
        double sumNorm = 0;
        for(double k = 0.01/fLc; k*fLc < 100; k*=aMultK){
        // for(double k = 1./fLc; k*fLc < 100; k*=aMultK){
            Ks.push_back(k);
            double normK = 0;
            // normK = k*k*k*pow(k*fLc,-gamma);

            if (k*fLc < 1)
                normK = k*k*k*pow(k*fLc,-gamma_large_scale);
            else
                normK = k*k*k*pow(k*fLc,-gamma);
            norms.push_back(normK);
            sumNorm += normK;
            fWaves.push_back(0);
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
        SmartPtr<MagneticField> mf = new TurbulentMF(r, lambdaMpc, Bgauss);
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
    void CellMF::GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        outValue[0] = 0;
        outValue[1] = fAmplitude*fCellDirection[int(fSize*x[2]/fCurrentMaxz)];
        outValue[2] = 0;
        // std::cout << fCurrentMaxz << '\t' << fLc << '\t' << x[2] << '\t' << int(fSize*x[2]/fCurrentMaxz) << '\t' << outValue[1] << '\n';
        // std::cin.get();
    }


    TableLOS_MF::TableLOS_MF(std::string filename)
    {
    	hid_t       file, dataset;         /* handles */
    	hid_t       datatype, dataspace;
    	hid_t       memspace;
    	H5T_class_t t_class;                 /* data type class */
    	H5T_order_t order;                 /* data order */
    	size_t      size;                  /*
					        * size of the data element
					        * stored in file
					        */
    	hsize_t     dimsm[1];              /* memory space dimensions */
    	hsize_t     dims_out[1];           /* dataset dimensions */
    	herr_t      status;
	
    	double         x[150000]; /* output buffer */
    	double         y[150000]; /* output buffer */
    	double         z[150000]; /* output buffer */
    	double 		   l_units[1];

    	double         Bx[150000]; /* output buffer */
    	double         By[150000]; /* output buffer */
    	double         Bz[150000]; /* output buffer */
    	double 		   B_units[1];

    	hsize_t      count[2];              /* size of the hyperslab in the file */
    	hsize_t      offset[2];             /* hyperslab offset in the file */
    	hsize_t      count_out[3];          /* size of the hyperslab in memory */
    	hsize_t      offset_out[3];         /* hyperslab offset in memory */
    	int          i, j, k, status_n, rank;

    	file = H5Fopen("B_8_7.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT);
    	// file = H5Fopen("Mrk501_h.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT);
    	dataset = H5Dopen2(file, "x_pos", H5P_DEFAULT);
   		dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, x);

    	dataset = H5Dopen2(file, "y_pos", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, y);

    	dataset = H5Dopen2(file, "z_pos", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, z);

    	dataset = H5Dopen2(file, "l_units", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, l_units);

    	dataset = H5Dopen2(file, "Bx", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, Bx);

    	dataset = H5Dopen2(file, "By", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, By);

    	dataset = H5Dopen2(file, "Bz", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, Bz);

    	dataset = H5Dopen2(file, "B_units", H5P_DEFAULT);
    	dataspace = H5Dget_space(dataset);    /* dataspace handle */
    	rank      = H5Sget_simple_extent_ndims(dataspace);
    	status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    	dimsm[0] = dims_out[0];
    	memspace = H5Screate_simple(1, dimsm, NULL);
    	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			     H5P_DEFAULT, B_units);


    	for (i=0; i<100; i++)
    		std::cout << x[i] << " " << y[i] << " " << z[i] << " " << l_units[0] << " "
    		          << Bx[i] << " " << By[i] << " " << Bz[i] << " " << B_units[0] << '\n';
    }

    void TableLOS_MF::GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        int index = 0;
        double fraction = 0;
        index = int(x[2]/MinVariabilityScale(aTime));
        fraction = x[2]/MinVariabilityScale(aTime) - index;
        outValue[0] = fBx[index] + fraction*(fBx[index + 1] - fBx[index]);
        outValue[1] = fBy[index] + fraction*(fBy[index + 1] - fBy[index]);
        outValue[2] = fBz[index] + fraction*(fBz[index + 1] - fBz[index]);
        // std::cout << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << MinVariabilityScale(aTime) << '\n';
        // std::cin.get();      
    }
    double TableLOS_MF::MinVariabilityScale(const CosmoTime &aTime) const
    {
        return 0.02*units.Mpc;
    }



    TableLOSTube_MF::TableLOSTube_MF()
    {
    	double tmp2;
        std::fstream fin("LOS_steps.dat");
        while (fin >> tmp2)
        {
            fL.push_back(tmp2*units.Mpc);
        }
        fin.close();

        std::vector<double> tmp;
        for(int j = 0; j < 5536; j++)
        	tmp.push_back(0);

        fin.open("test_Tube_z0_x.dat");
        for(int i = 0; i < 2500; i++)
        {
        	for(int j = 0; j < 5536; j++)
        	{
        		fin >> tmp[j];
        	}
        	fBx.push_back(tmp);
        }
        fin.close();

        fin.open("test_Tube_z0_y.dat");
        for(int i = 0; i < 2500; i++)
        {
        	for(int j = 0; j < 5536; j++)
        	{
        		fin >> tmp[j];
        	}
        	fBy.push_back(tmp);
        }
        fin.close();

        fin.open("test_Tube_z0_z.dat");
        for(int i = 0; i < 2500; i++)
        {
        	for(int j = 0; j < 5536; j++)
        	{
        		fin >> tmp[j];
        	}
        	fBz.push_back(tmp);
        }
        fin.close();
    }
    
    void TableLOSTube_MF::GetValueGauss(const double *x, const CosmoTime &aTime, std::vector<double> &outValue) const
    {
        int ix, iy, iz;
        double fx, fy, fz;
        double valx_1, valx_2, valx_3, valx_4, valx_5, valx_6, valx_7, valx_8; 
        double valxx_1, valxx_2, valxx_3, valxx_4;
        double valxxx_1, valxxx_2; 
        double valy_1, valy_2, valy_3, valy_4, valy_5, valy_6, valy_7, valy_8;
        double valyy_1, valyy_2, valyy_3, valyy_4;          
        double valyyy_1, valyyy_2; 
        double valz_1, valz_2, valz_3, valz_4, valz_5, valz_6, valz_7, valz_8;
        double valzz_1, valzz_2, valzz_3, valzz_4;  
        double valzzz_1, valzzz_2; 
        ix = int(x[0]/MinVariabilityScale(aTime));
        if (ix < 0)
        	ix--;
        iy = int(x[1]/MinVariabilityScale(aTime));
        if (iy < 0)
        	iy--;
        iz = int(x[2]/MinVariabilityScale(aTime));
        fx = x[0]/MinVariabilityScale(aTime) - ix;
        fy = x[1]/MinVariabilityScale(aTime) - iy;
        fz = x[2]/MinVariabilityScale(aTime) - iz;
        if (ix > 23)
        {
        	ix = 23;
        	fx = 0;
        }
        if (ix < -23)
        {
        	ix = -23;
        	fx = 0;
        }
        if (iy > 23)
        {
        	iy = 23;
        	fy = 0;
        }
        if (iy < -23)
        {
        	iy = -23;
        	fy = 0;
        }

        // x
        valx_1 = fBx[1200 + 49*iy - ix][iz];
        valx_2 = fBx[1200 + 49*iy - ix][iz + 1];
        valx_3 = fBx[1200 + 49*iy - (ix + 1)][iz];
        valx_4 = fBx[1200 + 49*iy - (ix + 1)][iz + 1];
        valx_5 = fBx[1200 + 49*(iy + 1) - (ix + 1)][iz];
        valx_6 = fBx[1200 + 49*(iy + 1) - (ix + 1)][iz + 1];
        valx_7 = fBx[1200 + 49*(iy + 1) - ix][iz];
        valx_8 = fBx[1200 + 49*(iy + 1) - ix][iz + 1];
        valxx_1 = valx_1 + fz*(valx_2 - valx_1);
        valxx_2 = valx_3 + fz*(valx_4 - valx_3);
        valxx_3 = valx_5 + fz*(valx_6 - valx_5);
        valxx_4 = valx_7 + fz*(valx_8 - valx_7);
        valxxx_1 = valxx_1 + fx*(valxx_2 - valxx_1);
        valxxx_2 = valxx_4 + fx*(valxx_3 - valxx_4);
        outValue[0] = valxxx_1 + fy*(valxxx_2 - valxxx_1);

        // y
        valy_1 = fBy[1200 + 49*iy - ix][iz];
        valy_2 = fBy[1200 + 49*iy - ix][iz + 1];
        valy_3 = fBy[1200 + 49*iy - (ix + 1)][iz];
        valy_4 = fBy[1200 + 49*iy - (ix + 1)][iz + 1];
        valy_5 = fBy[1200 + 49*(iy + 1) - (ix + 1)][iz];
        valy_6 = fBy[1200 + 49*(iy + 1) - (ix + 1)][iz + 1];
        valy_7 = fBy[1200 + 49*(iy + 1) - ix][iz];
        valy_8 = fBy[1200 + 49*(iy + 1) - ix][iz + 1];
        valyy_1 = valy_1 + fz*(valy_2 - valy_1);
        valyy_2 = valy_3 + fz*(valy_4 - valy_3);
        valyy_3 = valy_5 + fz*(valy_6 - valy_5);
        valyy_4 = valy_7 + fz*(valy_8 - valy_7);
        valyyy_1 = valyy_1 + fx*(valyy_2 - valyy_1);
        valyyy_2 = valyy_4 + fx*(valyy_3 - valyy_4);
        outValue[1] = valyyy_1 + fy*(valyyy_2 - valyyy_1);

        // z        
        valz_1 = fBz[1200 + 49*iy - ix][iz];
        valz_2 = fBz[1200 + 49*iy - ix][iz + 1];
        valz_3 = fBz[1200 + 49*iy - (ix + 1)][iz];
        valz_4 = fBz[1200 + 49*iy - (ix + 1)][iz + 1];
        valz_5 = fBz[1200 + 49*(iy + 1) - (ix + 1)][iz];
        valz_6 = fBz[1200 + 49*(iy + 1) - (ix + 1)][iz + 1];
        valz_7 = fBz[1200 + 49*(iy + 1) - ix][iz];
        valz_8 = fBz[1200 + 49*(iy + 1) - ix][iz + 1];
        valzz_1 = valz_1 + fz*(valz_2 - valz_1);
        valzz_2 = valz_3 + fz*(valz_4 - valz_3);
        valzz_3 = valz_5 + fz*(valz_6 - valz_5);
        valzz_4 = valz_7 + fz*(valz_8 - valz_7);
        valzzz_1 = valzz_1 + fx*(valzz_2 - valzz_1);
        valzzz_2 = valzz_4 + fx*(valzz_3 - valzz_4);
        outValue[2] = valzzz_1 + fy*(valzzz_2 - valzzz_1);     
    }
    double TableLOSTube_MF::MinVariabilityScale(const CosmoTime &aTime) const
    {
        return 0.02*units.Mpc;
    }
}