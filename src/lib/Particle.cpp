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
        "A2",
        "A3",
        "A4",
        "A5",
        "A6",
        "A7",
        "A8",
        "A9",
        "A10",
        "A11",
        "A12",
        "A13",
        "A14",
        "A15",
        "A16",
        "A17",
        "A18",
        "A19",
        "A20",
        "A21",
        "A22",
        "A23",
        "A24",
        "A25",
        "A26",
        "A27",
        "A28",
        "A29",
        "A30",
        "A31",
        "A32",
        "A33",
        "A34",
        "A35",
        "A36",
        "A37",
        "A38",
        "A39",
        "A40",
        "A41",
        "A42",
        "A43",
        "A44",
        "A45",
        "A46",
        "A47",
        "A48",
        "A49",
        "A50",
        "A51",
        "A52",
        "A53",
        "A54",
        "A55",
        "A56"
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
        logStr << ParticleNames[Type] << "(E=" << Energy/units.eV << "eV, " << Time.ToString() << ")";
        return logStr.str();
    }

    double Particle::Mass(ParticleType aType){
        if(aType<EndLightParticle)
            return MassesMeV[aType]/units.Eunit;
        if(aType<EndNuclei){
            double meanNucleonMassMeV = 0.5*(Particle::Mass(Proton)+Particle::Mass(Neutron));
            return meanNucleonMassMeV*Interactions::CNucleus::getA(aType);
        }
        NOT_IMPLEMENTED // treat other particles when added
        return 0.;
    }

    double Particle::ElectricCharge(ParticleType aType){
        if(aType<EndLightParticle)
            return ElectricCharges[aType];
        if(aType<EndNuclei)
            return ElectricCharges[Proton]*Interactions::CNucleus::getZ(aType);
        NOT_IMPLEMENTED // treat other particles when added
        return 0.;
    }

}
