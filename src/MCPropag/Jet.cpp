/*
 * Jet.cpp
 *
 *  Created on: Sept 15, 2021
 *      Author: dka
 */


#include <math.h>
#include <iostream>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

#include "Jet.h"

namespace MCPropag {

using namespace std;

bool JetStatic_isItFirstCall=true;

double coneZ_weight(double cos,double maxCos){
	if(cos>maxCos)return 1.;
	cerr<<"Warning: coneZ_weight returns zero weight"<<endl;
	return 0;
}
double GaussConeZ_weight(double setup_theta, double setup_phi, long double * pdir,double oneSigma){
	if(fabs(oneSigma)<1.e-14){
		if(JetStatic_isItFirstCall){
			JetStatic_isItFirstCall=false;
			cerr<<"Warning: You specified jet type as Gauss but one sigma specified as zero. Weight will be set as 1."<<endl;
			cerr<<"Warning: Maybe you want change --jet-param to non-zero value, change jet type to cone or von Miser-Fisher (vMF) type?"<<endl;
		}
		return 1.;
	}
	double Px = cos(setup_phi)*sin(setup_theta);
	double Py = sin(setup_phi)*sin(setup_theta);
	double Pz = cos(setup_theta);
	double scal=Px*pdir[0]+Py*pdir[1]+Pz*pdir[2];
	double res=gsl_ran_gaussian_pdf(fabs(1.-scal), oneSigma);
	if(res<0.)cerr<<"WTF! res="<<res<<endl;
	return res;
//	return exp(- (scal-1.)*(scal-1.)/(2.*oneSigma*oneSigma))/sqrt(6.283)/oneSigma;
}

double vMF_weight(double setup_theta, double setup_phi, long double * pdir, double kappa){
	if(fabs(kappa)<1.e-14)return 1.;
	double Px = cos(setup_phi)*sin(setup_theta);
	double Py = sin(setup_phi)*sin(setup_theta);
	double Pz = cos(setup_theta);
	double scal=Px*pdir[0]+Py*pdir[1]+Pz*pdir[2];
	return kappa*exp(kappa*(scal-1))/(2*M_PI*(1-exp(-2*kappa))); // Excuse me, wtf?
//	return (1-exp(-2*kappa)) * exp(-kappa*(scal-1))/(2.*kappa);
}


} /* namespace MCPropag */
