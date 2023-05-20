/*
 * Output.h
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


#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "Particle.h"
#include <string>
#include "Utils.h"
#include <fstream>

namespace mcray {

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
		RawOutput3D(std::string aDir, bool aOverwriteExisting = false, bool aSaveSourceZ=false, bool aSaveSourceE=false, bool aSaveId=false, bool aSaveZprod=false):
				RawOutput(aDir, false, aOverwriteExisting),fSaveSourceZ(aSaveSourceZ), fSaveSourceE(aSaveSourceE), fSaveId(aSaveId), fSaveZprod(aSaveZprod){ }
		void WriteHeader(std::ostream& aOut, const Particle& aFirstParticle);
		void Write(std::ostream& aOut, const Particle& aParticle);
	private:
		bool fSaveSourceZ;
		bool fSaveSourceE;
		bool fSaveId;
        bool fSaveZprod;
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

} /* namespace mcray */
#endif /* OUTPUT_H_ */
