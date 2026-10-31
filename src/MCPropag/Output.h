/*
 * Output.h
 *
 *  Created on: Oct 2, 2011
 *      Author: ok
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "Particle.h"
#include <string>
#include "Utils.h"
#include <fstream>

namespace MCPropag {

class IOutput : public ISmartReferencedObj{
public:
	virtual void AddParticle(const Particle& aParticle) = 0;
	virtual void Flush(bool aFinal) = 0;
};

class SpectrumOutput : public TSmartReferencedObj<IOutput> {
public:
	SpectrumOutput(std::string aFile, double aEmin, double aStep);
	virtual ~SpectrumOutput();
	void AddParticle(const Particle& aParticle);
	void Flush(bool aFinal);
private:
	Vector			fSpectra[ParticleTypeEOF];
	std::vector<int> fNumberOfParticles[ParticleTypeEOF];
	double			fEmin;
	double			fLogStep;
	std::string 	fFile;
};

class RawOutput : public TSmartReferencedObj<IOutput> {
public:
	RawOutput(std::string aDir, bool aCheckRepetitions, bool aOverwriteExisting = false);
	virtual ~RawOutput();
	void AddParticle(const Particle& aParticle);
	void Flush(bool aFinal);
	void SetOutputDir(std::string aDir, bool aOverwriteExisting = false);
	virtual void WriteHeader(std::ostream& aOut, const Particle& aFirstParticle);
	virtual void Write(std::ostream& aOut, const Particle& aParticle);

	//TODO: delete this method when intermediate output is implemented via multiple IResult instances in PropagEngine
	void CopyToStack(class ParticleStack& aStack);
private:
	void LookForRepetitions(std::vector<Particle>* aParticles);
	std::vector<Particle>*	fParticles;
    std::vector<bool> 	    fHeaderWritten;
	std::string				fDir;
	bool					fCheckRepetitions;
    bool                    fFlushInProgress;
};

	class RawOutput3D : public RawOutput{
	public:
		RawOutput3D(std::string aDir, bool aOverwriteExisting = false, bool aSaveSourceZ=false, bool aSaveSourceE=false, bool aSaveId=false, bool aTargetLogger=false,bool aSaveSourceDir=false):
				RawOutput(aDir, false, aOverwriteExisting),fSaveSourceZ(aSaveSourceZ), fSaveSourceE(aSaveSourceE), fSaveId(aSaveId),
				fTargetLogger(aTargetLogger),fSaveSourceDir(aSaveSourceDir){ }
		void WriteHeader(std::ostream& aOut, const Particle& aFirstParticle);
		void Write(std::ostream& aOut, const Particle& aParticle);
	private:
		bool fSaveSourceZ;
		bool fSaveSourceE;
		bool fSaveId;
		bool fTargetLogger;
		bool fSaveSourceDir;
		double normT;
	};

	class TotalEnergyOutput : public TSmartReferencedObj<IOutput> {
	public:
		TotalEnergyOutput(){memset(fE,0,ParticleTypeEOF*sizeof(double));fTotE=0.;}
		void AddParticle(const Particle& aParticle) {
			fTotE += aParticle.Energy*aParticle.Weight;
			fE[aParticle.Type] += aParticle.Energy*aParticle.Weight;
		}
		void Flush(bool aFinal){}
		double TotalE() const { return fTotE; }
		double TotalE(ParticleType aType) const { return fE[aType]; }
	private:
		double	fTotE;
		double	fE[ParticleTypeEOF];
	};

} /* namespace MCPropag */
#endif /* OUTPUT_H_ */
