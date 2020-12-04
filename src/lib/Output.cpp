/*
 * Output.cpp
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


#include "Output.h"
#include<sys/stat.h>
#include <map>
#include <utility>
#include "ParticleStack.h"

namespace mcray {

SpectrumOutput::SpectrumOutput(std::string aFile, double aEmin, double aStep):
fEmin(aEmin),
fLogStep(log(aStep)),
fFile(aFile)
{
	fEmin /= sqrt(aStep);//center bins in Emin*aStep^N
	//std::ofstream	fileOut;
	//fileOut.open(fFile.c_str());//make sure file can be created and delete the contents of old file if it exists
	//if(!fileOut)
		//Exception::Throw("Failed to open " + fFile + " for writing");
	//fileOut.close();
}

SpectrumOutput::~SpectrumOutput() {

}

void SpectrumOutput::AddParticle(const Particle& aParticle)
{
	int iBin = log(aParticle.Energy/fEmin)/fLogStep;
	if(iBin<0)
		return;//skip output of particles with E<Emin
	Vector& spec = fSpectra[aParticle.Type];
	std::vector<int>& numbers = fNumberOfParticles[aParticle.Type];
	while((int)spec.size()<=iBin)
	{
		spec.push_back(0);
		numbers.push_back(0);
	}
	spec[iBin] += aParticle.Weight;
	numbers[iBin] += 1;
}

void SpectrumOutput::Flush(bool aFinal)
{
    if(!aFinal)
        return;//intermedeate output is not supported
	std::ofstream	fileOut;
	fileOut.open(fFile.c_str());//replace previous content
	//fileOut.open(fFile.c_str(), std::ios::app);//don't remove previous contents
	if(!fileOut)
		Exception::Throw("Failed to open " + fFile + " for writing");
	//fileOut.seekp(std::ios::beg);//every next spectrum output is saved in the beginning of the file
	unsigned int nBins = 0;
	for(int p=0; p<ParticleTypeEOF; p++)
	{
		unsigned int size = fSpectra[p].size();
		if(size>nBins)
			nBins = size;
	}
	double mult = exp(fLogStep);
	double E=fEmin*sqrt(mult);
	double deltaEmult = sqrt(mult) - 1./sqrt(mult);
	for(unsigned int iE=0; iE<nBins; iE++, E*=mult)
	{
		fileOut << E*units.Eunit*1e6;//eV
		for(unsigned int p=0; p<ParticleTypeEOF; p++)
		{
			double intFlux = 0;
			int	nParticles = 0;
			if(fSpectra[p].size()>iE)
			{
				intFlux = fSpectra[p][iE];
				nParticles = fNumberOfParticles[p][iE];
			}
			fileOut << "\t" << intFlux*deltaEmult*E << "\t" << nParticles;//print dN(E)/dE * E^2
		}
		fileOut << '\n';
	}
	fileOut << "\n\n";/// delimiters between consequent spectra outputs
	fileOut.close();
}

RawOutput::RawOutput(std::string aDir, bool aCheckRepetitions, bool aOverwriteExisting):
		fParticles(0),
		fDir(aDir),
		fCheckRepetitions(aCheckRepetitions),
        fHeaderWritten(ParticleTypeEOF,false),
        fFlushInProgress(false)
{
	SetOutputDir(aDir, aOverwriteExisting);
	fParticles = new std::vector<Particle>[ParticleTypeEOF];
}

RawOutput::~RawOutput()
{
	delete[] fParticles;
}

void RawOutput::SetOutputDir(std::string aDir, bool aOverwriteExisting)
{
	fDir = aDir;
	if(mkdir(fDir.c_str(),0777))
	{
		if(aOverwriteExisting)
		{
			system(("rm -r " + aDir).c_str());
			if(!mkdir(fDir.c_str(),0777))
				return;
		}
		Exception::Throw("Failed to create directory " + aDir);
	}
}

void RawOutput::CopyToStack(ParticleStack& aStack)
{
	for(int p=0; p<ParticleTypeEOF; p++)
	{
		std::vector<Particle>& particles = fParticles[p];
		if(particles.size())
			aStack.AddSecondaries(particles);
	}
}

void RawOutput::AddParticle(const Particle& aParticle)
{
	fParticles[aParticle.Type].push_back(aParticle);
}

void RawOutput::Flush(bool aFinal)
{
	if(fFlushInProgress)
		Exception::Throw("RawOutput::Flush(): was invoked from more than one thread");
    fFlushInProgress = true;
    if(fCheckRepetitions && aFinal)
        LookForRepetitions(fParticles);


	for(int p=0; p<ParticleTypeEOF; p++)
	{
		std::vector<Particle>& particles = fParticles[p];
		if(particles.size()==0)
			continue;
		std::string fileName = fDir + "/" + Particle::Name((ParticleType)p);
		std::ofstream file;
		file.open(fileName.c_str(), std::ios::app);
		if(!file) {
            fFlushInProgress = false;
            Exception::Throw("Failed to open " + fileName + " for writing");
        }
		std::vector<Particle>::const_iterator pit = particles.begin();
        if(!fHeaderWritten[p]) {
            WriteHeader(file, *pit);
            file << "\n";
            fHeaderWritten[p] = true;
        }
		for(; pit != particles.end(); pit++)
		{
			Write(file, *pit);
			file << "\n";
		}
		file.close();
		particles.clear();
	}
    fFlushInProgress = false;
}

void RawOutput::WriteHeader(std::ostream& aOut, const Particle& aFirstParticle){
	aOut << "# E/eV\tweight\tNinteractions\tz_fin\tsinBeta\tD_lastDeflection/Mpc\tD_source\tz_source";
}

void RawOutput::Write(std::ostream& aOut, const Particle& aParticle){
	double beta = aParticle.DeflectionAngle();
	double sinBeta = (beta<0.5*M_PI)?sin(aParticle.DeflectionAngle()):1.;//to enable correct filtering based on sin(beta) set sin=1 for particles deflected by angles larger than Pi/2
	double d = cosmology.z2d(aParticle.SourceParticle->Time.z())/units.Mpc;
	double l = cosmology.z2d(aParticle.LastDeflectionTime().z())/units.Mpc;
	aOut << aParticle.Energy*units.Eunit*1e6 << "\t" << aParticle.Weight << "\t"
	<< aParticle.Ninteractions << "\t" << aParticle.Time.z() <<  "\t" //used for debugging
	<< sinBeta << "\t" << l << "\t" << d << "\t" << aParticle.SourceParticle->Time.z(); //used to filter on base of PSF and jet angle or on base of source Z
}

void RawOutput::LookForRepetitions(std::vector<Particle>* aParticles)
{
	std::map<float, int> energies;
	for(int p=0; p<ParticleTypeEOF; p++)
	{
		std::vector<Particle>& particles = aParticles[p];
		for(std::vector<Particle>::const_iterator pit = particles.begin(); pit != particles.end(); pit++)
		{
			if(!pit->Ninteractions)//skipping particles which did not interact
				continue;
			float log10E = (float)log10(pit->Energy*units.Eunit*1e6);
			energies[log10E]++;
		}
	}

	std::ofstream file;
	for (std::map<float,int>::iterator it = energies.begin(); it != energies.end(); it++)
	{
		if(it->second > 1)
		{
			if(!file.is_open())
			{
				std::cerr << "Repetitions found" << std::endl;
				std::string fileName = fDir + "/repetitions";
				file.open(fileName.c_str(),std::ios::out);
				if(!file.is_open())
					Exception::Throw("Failed to create " + fileName);
			}
			file << it->first << "\t" << it->second << "\n";
		}
	}
	if(file.is_open())
		file.close();
}

	void RawOutput3D::Write(std::ostream& aOut, const Particle& aParticle){
		double T = aParticle.Time.t()-aParticle.SourceParticle->Time.t();
		aOut << aParticle.Energy/units.eV << "\t" << aParticle.Weight << "\t"
		<< aParticle.X[0]/T << "\t" << aParticle.X[1]/T << "\t" << 1.-aParticle.X[2]/T << "\t"
		<< aParticle.Pdir[0] << "\t" << aParticle.Pdir[1] << "\t" << 1.-aParticle.Pdir[2] << "\t"
		<< aParticle.fCascadeProductionTime.z() << "\t" << aParticle.Ninteractions;
		if(fSaveSourceZ)
			aOut << "\t" << aParticle.SourceParticle->Time.z();
		if(fSaveSourceE)
			aOut << "\t" << aParticle.SourceParticle->Energy/units.eV;
		if(fSaveId){
			aOut << "\t" << aParticle.id;
		}
	}

	void RawOutput3D::WriteHeader(std::ostream& aOut, const Particle& aFirstParticle){
		aOut << "# " << Particle::Name(aFirstParticle.Type);
		if(!fSaveSourceZ){
			aOut << "\tz_src = "<< aFirstParticle.SourceParticle->Time.z();
		}
		if(!fSaveSourceE){
			aOut << "\tE_src/eV = "<< aFirstParticle.SourceParticle->Energy/units.eV;
		}
		aOut << "\t" << "T/Mpc = " << (aFirstParticle.Time.t()-aFirstParticle.SourceParticle->Time.t())/units.Mpc << "\n#\n";
		aOut << "# E/eV\tweight\tx/T\ty/T\t1-(z/T)\tPx/P\tPy/P\t1-(Pz/P)\tz_cascade_production\tN_interactions";
		if(fSaveSourceZ){
			aOut << "\tz_src";
		}
		if(fSaveSourceE){
			aOut << "\tE_src/eV";
		}
		if(fSaveId){
			aOut << "\tId";
		}
	}

} /* namespace mcray */
