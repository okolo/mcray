/*
 * Jet.cpp
 *
 *  Created on: Sept 15, 2021
 *      Author: dka
 */


#include <math.h>
#include <iostream>

#include "Jet.h"

namespace MCPropag {

using namespace std;

double coneZ_weight(double cos,double maxCos){
	if(cos>maxCos)
		return 1./maxCos;
	cerr<<"Warning: coneZ_weight returns zero weight"<<endl;
	return 0;
}
double GaussConeZ_weight(double cos,double oneSigma){
	return 1./sqrt(2.*M_PI)/oneSigma*exp(- cos*cos/(2.*oneSigma*oneSigma));
}

} /* namespace MCPropag */
