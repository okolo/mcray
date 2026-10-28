/*
 * Particle.cpp
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


#include "Particle.h"
#include "Nucleus.h"
#include <iomanip>

namespace mcray
{

const double Particle::MassesMeV[EndLightParticle] =
{// order of masses in this array must follow the order of particles in the enumeration ParticleType
//#ifdef ELMAG_TEST
		0.511, //Electron,
		0.511,//Positron,
//#else
		//0.510998928, //Electron,
		//0.510998928,//Positron,
//#endif
        0,//Photon,
        0,//NeutrinoE,
        0,//NeutrinoM,
        0,//NeutrinoT,
        0,//NeutrinoAE,
        0,//NeutrinoAM,
        0,//NeutrinoAT,
        939.565378,//Neutron,
        938.272046//Proton,

};

const int Particle::ElectricCharges[EndLightParticle] =
{// order of masses in this array must follow the order of particles in the enumeration ParticleType
		-1, //Electron,
		1,//Positron,
        0,//Photon,
        0,//NeutrinoE,
        0,//NeutrinoM,
        0,//NeutrinoT,
        0,//NeutrinoAE,
        0,//NeutrinoAM,
        0,//NeutrinoAT,
        0,//Neutron,
        1//Proton,
};

const char* Particle::ParticleNames[ParticleTypeEOF] =
{
	    "electron",
	    "positron",
	    "photon",
	    "neutrinoE",
	    "neutrinoM",
	    "neutrinoT",
	    "neutrinoAE",
	    "neutrinoAM",
	    "neutrinoAT",
	    "neutron",
	    "proton",
        "nucleus"
};

    int Particle::fNextCustomDataIndex = 0;
    unsigned long long Particle::fMaxId = 0;

    void Particle::PropagateFreely(cosmo_time dt)
    {
        Time += dt;
        coord_type b=beta();
        for(int i=0; i<3; i++)
            X[i]+=(b*Pdir[i]*dt);
    }

    void Particle::GenerateId(){
#pragma omp critical (Particle)
        id = ++fMaxId;
    }

    int Particle::ReserveInteractionDataSlot()//returns index to access and store interaction data or -1 if no space left
    {
        int result = -1;
#pragma omp critical (Particle)
        {
            if(fNextCustomDataIndex < ParticleCustomDataSize)
                result = fNextCustomDataIndex;
            fNextCustomDataIndex++;
        }
        return result;
    }

    std::string Particle::ToString(){
        std::ostringstream logStr;
        logStr << Name() << "(E=" << Energy/units.eV << "eV, " << Time.ToString() << ")";
        return logStr.str();
    }

    void Particle::ValidateIdentity() const {
        if (Type < Electron || Type >= ParticleTypeEOF)
            Exception::Throw("Invalid particle type");
        if (Type == Nucleus) {
            if (fA < 2 || !fZ || fZ > fA || !std::isfinite(fNuclearMass) ||
                fNuclearMass <= 0. || fMassModel == NuclearMassModel::Unset)
                Exception::Throw("Invalid or unset nuclear A,Z");
        } else if (fA || fZ || fNuclearMass || fMassModel != NuclearMassModel::Unset) Exception::Throw("Particle type changed without clearing nuclear identity; use SetType");
    }
    void Particle::SetType(ParticleType type) {
        if (type < Electron || type >= Nucleus)
            Exception::Throw("Use a nucleus setter to establish A,Z");
        Type = type; fA = fZ = 0; fNuclearMass = 0.; fMassModel = NuclearMassModel::Unset;
    }
    void Particle::SetPsbNucleus(unsigned int a) {
        if (a < 2 || a > Interactions::CNucleus::MaxA)
            Exception::Throw("Unsupported PSB mass number");
        SetNucleus(a, Interactions::CNucleus::getZ(a),
                   a * 0.5 * (Mass(Proton) + Mass(Neutron)), NuclearMassModel::PSB);
    }
    void Particle::SetNucleus(unsigned int a, unsigned int z, double mass, NuclearMassModel model) {
        if (a < 2 || !z || z > a || !std::isfinite(mass) || mass <= 0. ||
            (model != NuclearMassModel::PSB && model != NuclearMassModel::TALYS))
            Exception::Throw("Invalid nuclear identity or mass");
        Type = Nucleus; fA = a; fZ = z; fNuclearMass = mass; fMassModel = model;
    }
    unsigned int Particle::MassNumber() const {
        ValidateIdentity();
        return Type == Nucleus ? fA : (Type == Proton || Type == Neutron ? 1 : 0);
    }
    unsigned int Particle::AtomicNumber() const {
        ValidateIdentity();
        return Type == Nucleus ? fZ : (Type == Proton ? 1 : 0);
    }
    double Particle::Mass() const {
        ValidateIdentity();
        return Type == Nucleus ? fNuclearMass : Mass(Type);
    }
    int Particle::ElectricCharge() const {
        ValidateIdentity();
        return Type == Nucleus ? fZ : ElectricCharge(Type);
    }
    std::string Particle::Name() const {
        ValidateIdentity();
        if (Type != Nucleus) return Name(Type);
        std::ostringstream name;
        name << 'A' << std::setfill('0') << std::setw(2) << fA
             << 'Z' << std::setw(2) << fZ;
        return name.str();
    }

    double Particle::Mass(ParticleType aType){
        if(aType>=Electron && aType<EndLightParticle)
            return MassesMeV[aType]/units.Eunit;
        Exception::Throw("Nuclear mass/charge requires a Particle instance with A,Z");
        return 0.;
    }

    double Particle::ElectricCharge(ParticleType aType){
        if(aType>=Electron && aType<EndLightParticle)
            return ElectricCharges[aType];
        Exception::Throw("Nuclear mass/charge requires a Particle instance with A,Z");
        return 0.;
    }

}
