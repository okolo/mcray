/*
 * Observer.h
 *
 *  Created on: Oct 4, 2011
 *      Author: ok
 */

#ifndef OBSERVER_H_
#define OBSERVER_H_

#include "Particle.h"
#include "Cosmology.h"

namespace MCPropag {

class Observer{
public:
//	Observer();
	Observer(int type=-1, double x=0,double y=0,double z=0, double R=0);
	bool weInTarget(const Particle& aParticle) const;
	cosmo_time patchDtIfNeed(const Particle& aParticle,cosmo_time maxDt);
	int type;
private:
	double x,y,z,R,R2,R3,Rdt;
};

} /* namespace MCPropag */
#endif /* OBSERVER_H_ */
