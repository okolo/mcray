/*
 * Jet.h
 *
 *  Created on: Sept 15, 2021
 *      Author: dka
 */

#ifndef JET_H_
#define JET_H_

namespace MCPropag {

double coneZ_weight(double cos,double maxCos);
double GaussConeZ_weight(double setup_theta, double setup_phi, long double * pdir,double oneSigma);
double vMF_weight(double setup_theta, double setup_phi, long double * pdir, double kappa);

} /* namespace MCPropag */
#endif /* JET_H_ */
