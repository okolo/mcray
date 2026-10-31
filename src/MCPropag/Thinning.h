/* 
 * File:   Thinning.h
 * Author: ok
 *
 * Created on August 26, 2011, 7:01 PM
 */

#ifndef THINNING_H
#define	THINNING_H

#include <vector>
#include "Particle.h"
#include "Randomizer.h"

namespace MCPropag
{

class IThinning {
public:
	//must be thread-safe
    virtual void Run(std::vector<Particle>& aParticles, Randomizer& aRandomizer) = 0;
};

	enum FilterMode{
		BlackList,
		WhiteList
	};

class EnergyBasedThinning : public  IThinning {
public:
	EnergyBasedThinning(double aAlpha, FilterMode aFilterMode=WhiteList);
    void Run(std::vector<Particle>& aParticles, Randomizer& aRandomizer);
    virtual ~EnergyBasedThinning(){}
	void EnableFor(ParticleType aType, bool aDoEnable = true) { fEnableFor[aType] = aDoEnable;}
private:
    double				fAlpha;
	std::vector<bool>	fEnableFor;
};


// Added by AK to treat low energy photons from EM cascades
class PhotonThinning : public  IThinning {
public:
	// thinning works for energies below fEnergy
	PhotonThinning(double aAlpha, double aThreshold);
    void Run(std::vector<Particle>& aParticles, Randomizer& aRandomizer);
    virtual ~PhotonThinning(){}
private:
    double				fAlpha;
    double				fEnergy;  // added by AK
};

}
#endif	/* THINNING_H */

