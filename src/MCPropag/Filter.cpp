#include "Filter.h"

namespace MCPropag
{

	JetOutputFilter::JetOutputFilter(double aEmin, double aJetOpenningAngle, double aObservationAngle) :
			fEmin(aEmin),
			fJetOpenningAngle(aJetOpenningAngle),
			fObservationAngle(aObservationAngle)
	{
	}

	bool JetOutputFilter::Pass(const Particle& aParticle) const
	{
		if(aParticle.Energy<=fEmin)
			return false;
		if(fObservationAngle<M_PI_2 && aParticle.ObservationAngle() > fObservationAngle)
			return false;
		if(fJetOpenningAngle<M_PI && aParticle.JetOpenningAngle() > fJetOpenningAngle)
			return false;
		return true;
	}

	EnergyBasedFilter::EnergyBasedFilter(double aEmin, double aMaxDeflectionAngle) :
			fEmin(aEmin),
			fMaxDeflectionAngle(aMaxDeflectionAngle)
	{
	}

	bool EnergyBasedFilter::Pass(const Particle& aParticle) const
	{
		if (aParticle.Energy < fEmin)
		{
			return false;
		}	
		if(fMaxDeflectionAngle<M_PI && aParticle.DeflectionAngle() > fMaxDeflectionAngle)
		{
			return false;
		}
		return true;
	}

	TypeSpecificEnergyBasedFilter::TypeSpecificEnergyBasedFilter(std::vector<int> aFilterParticleType, 
		std::vector<double> aFilterParticleEnergy, double aMaxDeflectionAngle) :
			fMaxDeflectionAngle(aMaxDeflectionAngle)
	{
		for (int i = 0; i < aFilterParticleType.size(); i++)
		{
			fFilterParticleType.push_back(aFilterParticleType[i]);
			fFilterParticleEnergy.push_back(aFilterParticleEnergy[i]);
		}
	}

	bool TypeSpecificEnergyBasedFilter::Pass(const Particle& aParticle) const
	{
		for (int i = 0; i < fFilterParticleType.size(); i++)
		{
			if (aParticle.Type == fFilterParticleType[i])
				if (aParticle.Energy < fFilterParticleEnergy[i])
					return false;
		}		
		if(fMaxDeflectionAngle<M_PI && aParticle.DeflectionAngle() > fMaxDeflectionAngle)
		{
			return false;
		}
		return true;
	}

	JetFilter::JetFilter(double aEmin, double aJetOpenningAngle, double aObservationAngle) :
			fEmin(aEmin),
			fJetOpenningAngleSin(aJetOpenningAngle<M_PI_2 ? sin(aJetOpenningAngle) : 1.),
			fObservationAngleCos(cos(aObservationAngle)),
			fObservationAngleSin(sin(aObservationAngle))
	{
	}

	bool JetFilter::Pass(const Particle& aParticle) const
	{
		if(aParticle.Energy<=fEmin)
			return false;
		if(aParticle.DeflectionAngle()>0.5*M_PI && fJetOpenningAngleSin<1.)
			return false;//exclude particles deflected by more than Pi/2
		double sinBeta = sin(aParticle.DeflectionAngle());
		double d = cosmology.z2d(aParticle.SourceParticle->Time.z());
		double l = cosmology.z2d(aParticle.LastDeflectionTime().z());

		if(sinBeta > d/l*fJetOpenningAngleSin)
			return false;//particle can not reach observer from any angle

		double r = sqrt(d*d+l*l-2.*d*l*fObservationAngleCos);//distance(at z=0) from the point of last deflection time to the source
															//assuming that particle arrives at the edge of PSF
		if(sinBeta > d/r*fObservationAngleSin)
			return false;//particle will reach observer from outside PSF

		return true;
	}

	ParticleTypeFilter::ParticleTypeFilter(const ParticleType aParticlesToDiscard[]) {
		for(int i=0; i<ParticleTypeEOF; i++){
			fFilter[i] = true;
		}
		if(aParticlesToDiscard)
			for(int i=0; aParticlesToDiscard[i]!=ParticleTypeEOF; i++) {
				fFilter[aParticlesToDiscard[i]] = false;
			}
	}
}
