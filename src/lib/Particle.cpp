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
#include <tuple>

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

const int Particle::NucleonStructure[EndRealNuclei-StartRealNuclei][2] =
{// order of masses in this array must follow the order of particles in the enumeration ParticleType
        //{Z, N}
        {1, 1},//"H2",
        {1, 2},//"H3",
        {2, 1},//"He3",
        {2, 2},//"He4",
        {3, 3},//"Li6",
        {3, 4},//"Li7",
        {4, 3},//"Be7",
        {4, 5},//"Be9",
        {4, 6},//"Be10",
        {4, 7},//"Be11",
        {5, 5},//"B10",
        {5, 6},//"B11",
        {6, 4},//"C10",
        {6, 5},//"C11",
        {6, 6},//"C12",
        {6, 7},//"C13",
        {6, 8},//"C14",
        {6, 9},//"C15",
        {6, 10},//"C16",
        {7, 6},//"N13",
        {7, 7},//"N14",
        {7, 8},//"N15",
        {7, 9},//"N16",
        {7, 10},//"N17",
        {8, 6},//"O14",
        {8, 7},//"O15",
        {8, 8},//"O16",
        {8, 9},//"O17",
        {8, 10},//"O18",
        {8, 11},//"O19",
        {8, 12},//"O20",
        {8, 13},//"O21",
        {8, 14},//"O22",
        {9, 8},//"F17",
        {9, 9},//"F18",
        {9, 10},//"F19",
        {9, 11},//"F20",
        {9, 12},//"F21",
        {9, 13},//"F22",
        {9, 14},//"F23",
        {10, 8},//"Ne18",
        {10, 9},//"Ne19",
        {10, 10},//"Ne20",
        {10, 11},//"Ne21",
        {10, 12},//"Ne22",
        {10, 13},//"Ne23",
        {10, 14},//"Ne24",
        {11, 10},//"Na21",
        {11, 11},//"Na22",
        {11, 12},//"Na23",
        {11, 13},//"Na24",
        {11, 14},//"Na25",
        {12, 10},//"Mg22",
        {12, 11},//"Mg23",
        {12, 12},//"Mg24",
        {12, 13},//"Mg25",
        {12, 14},//"Mg26",
        {12, 15},//"Mg27",
        {12, 16},//"Mg28",
        {13, 11},//"Al24",
        {13, 12},//"Al25",
        {13, 13},//"Al26",
        {13, 14},//"Al27",
        {13, 15},//"Al28",
        {13, 16},//"Al29",
        {13, 17},//"Al30",
        {14, 12},//"Si26",
        {14, 13},//"Si27",
        {14, 14},//"Si28",
        {14, 15},//"Si29",
        {14, 16},//"Si30",
        {14, 17},//"Si31",
        {14, 18},//"Si32",
        {14, 19},//"Si33",
        {14, 20},//"Si34",
        {15, 14},//"P29",
        {15, 15},//"P30",
        {15, 16},//"P31",
        {15, 17},//"P32",
        {15, 18},//"P33",
        {15, 19},//"P34",
        {15, 20},//"P35",
        {15, 21},//"P36",
        {15, 22},//"P37",
        {16, 15},//"S31",
        {16, 16},//"S32",
        {16, 17},//"S33",
        {16, 18},//"S34",
        {16, 19},//"S35",
        {16, 20},//"S36",
        {16, 21},//"S37",
        {16, 22},//"S38",
        {16, 23},//"S39",
        {16, 24},//"S40",
        {16, 25},//"S41",
        {17, 16},//"Cl33",
        {17, 17},//"Cl34",
        {17, 18},//"Cl35",
        {17, 19},//"Cl36",
        {17, 20},//"Cl37",
        {17, 21},//"Cl38",
        {17, 22},//"Cl39",
        {17, 23},//"Cl40",
        {17, 24},//"Cl41",
        {17, 25},//"Cl42",
        {17, 26},//"Cl43",
        {18, 17},//"Ar35",
        {18, 18},//"Ar36",
        {18, 19},//"Ar37",
        {18, 20},//"Ar38",
        {18, 21},//"Ar39",
        {18, 22},//"Ar40",
        {18, 23},//"Ar41",
        {18, 24},//"Ar42",
        {18, 25},//"Ar43",
        {18, 26},//"Ar44",
        {18, 27},//"Ar45",
        {18, 28},//"Ar46",
        {19, 19},//"K38",
        {19, 20},//"K39",
        {19, 21},//"K40",
        {19, 22},//"K41",
        {19, 23},//"K42",
        {19, 24},//"K43",
        {19, 25},//"K44",
        {19, 26},//"K45",
        {19, 27},//"K46",
        {19, 28},//"K47",
        {19, 29},//"K48",
        {19, 30},//"K49",
        {20, 20},//"Ca40",
        {20, 21},//"Ca41",
        {20, 22},//"Ca42",
        {20, 23},//"Ca43",
        {20, 24},//"Ca44",
        {20, 25},//"Ca45",
        {20, 26},//"Ca46",
        {20, 27},//"Ca47",
        {20, 28},//"Ca48",
        {20, 29},//"Ca49",
        {20, 30},//"Ca50",
        {21, 22},//"Sc43",
        {21, 23},//"Sc44",
        {21, 24},//"Sc45",
        {21, 25},//"Sc46",
        {21, 26},//"Sc47",
        {21, 27},//"Sc48",
        {21, 28},//"Sc49",
        {21, 29},//"Sc50",
        {21, 30},//"Sc51",
        {22, 22},//"Ti44",
        {22, 23},//"Ti45",
        {22, 24},//"Ti46",
        {22, 25},//"Ti47",
        {22, 26},//"Ti48",
        {22, 27},//"Ti49",
        {22, 28},//"Ti50",
        {22, 29},//"Ti51",
        {22, 30},//"Ti52",
        {23, 24},//"V47",
        {23, 25},//"V48",
        {23, 26},//"V49",
        {23, 27},//"V50",
        {23, 28},//"V51",
        {23, 29},//"V52",
        {23, 30},//"V53",
        {24, 24},//"Cr48",
        {24, 25},//"Cr49",
        {24, 26},//"Cr50",
        {24, 27},//"Cr51",
        {24, 28},//"Cr52",
        {24, 29},//"Cr53",
        {24, 30},//"Cr54",
        {25, 26},//"Mn51",
        {25, 27},//"Mn52",
        {25, 28},//"Mn53",
        {25, 29},//"Mn54",
        {25, 30},//"Mn55",
        {26, 26},//"Fe52",
        {26, 27},//"Fe53",
        {26, 28},//"Fe54",
        {26, 29},//"Fe55",
        {26, 30},//"Fe56",
        {26, 31},//"Fe57",
        {26, 32},//"Fe58",
        {26, 33},//"Fe59",
        {26, 34},//"Fe60",
        {26, 35}//"Fe61"
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
        "A56",
        "H2",
        "H3",
        "He3",
        "He4",
        "Li6",
        "Li7",
        "Be7",
        "Be9",
        "Be10",
        "Be11",
        "B10",
        "B11",
        "C10",
        "C11",
        "C12",
        "C13",
        "C14",
        "C15",
        "C16",
        "N13",
        "N14",
        "N15",
        "N16",
        "N17",
        "O14",
        "O15",
        "O16",
        "O17",
        "O18",
        "O19",
        "O20",
        "O21",
        "O22",
        "F17",
        "F18",
        "F19",
        "F20",
        "F21",
        "F22",
        "F23",
        "Ne18",
        "Ne19",
        "Ne20",
        "Ne21",
        "Ne22",
        "Ne23",
        "Ne24",
        "Na21",
        "Na22",
        "Na23",
        "Na24",
        "Na25",
        "Mg22",
        "Mg23",
        "Mg24",
        "Mg25",
        "Mg26",
        "Mg27",
        "Mg28",
        "Al24",
        "Al25",
        "Al26",
        "Al27",
        "Al28",
        "Al29",
        "Al30",
        "Si26",
        "Si27",
        "Si28",
        "Si29",
        "Si30",
        "Si31",
        "Si32",
        "Si33",
        "Si34",
        "P29",
        "P30",
        "P31",
        "P32",
        "P33",
        "P34",
        "P35",
        "P36",
        "P37",
        "S31",
        "S32",
        "S33",
        "S34",
        "S35",
        "S36",
        "S37",
        "S38",
        "S39",
        "S40",
        "S41",
        "Cl33",
        "Cl34",
        "Cl35",
        "Cl36",
        "Cl37",
        "Cl38",
        "Cl39",
        "Cl40",
        "Cl41",
        "Cl42",
        "Cl43",
        "Ar35",
        "Ar36",
        "Ar37",
        "Ar38",
        "Ar39",
        "Ar40",
        "Ar41",
        "Ar42",
        "Ar43",
        "Ar44",
        "Ar45",
        "Ar46",
        "K38",
        "K39",
        "K40",
        "K41",
        "K42",
        "K43",
        "K44",
        "K45",
        "K46",
        "K47",
        "K48",
        "K49",
        "Ca40",
        "Ca41",
        "Ca42",
        "Ca43",
        "Ca44",
        "Ca45",
        "Ca46",
        "Ca47",
        "Ca48",
        "Ca49",
        "Ca50",
        "Sc43",
        "Sc44",
        "Sc45",
        "Sc46",
        "Sc47",
        "Sc48",
        "Sc49",
        "Sc50",
        "Sc51",
        "Ti44",
        "Ti45",
        "Ti46",
        "Ti47",
        "Ti48",
        "Ti49",
        "Ti50",
        "Ti51",
        "Ti52",
        "V47",
        "V48",
        "V49",
        "V50",
        "V51",
        "V52",
        "V53",
        "Cr48",
        "Cr49",
        "Cr50",
        "Cr51",
        "Cr52",
        "Cr53",
        "Cr54",
        "Mn51",
        "Mn52",
        "Mn53",
        "Mn54",
        "Mn55",
        "Fe52",
        "Fe53",
        "Fe54",
        "Fe55",
        "Fe56",
        "Fe57",
        "Fe58",
        "Fe59",
        "Fe60",
        "Fe61"
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
        if(aType<EndRealNuclei)
            return NucleonStructure[aType-StartRealNuclei][0]*Particle::Mass(Proton)+NucleonStructure[aType-StartRealNuclei][1]*Particle::Mass(Neutron);
        NOT_IMPLEMENTED // treat other particles when added
        return 0.;
    }

        double Particle::AtomicMass(ParticleType aType){
        if(aType<=EndRealNuclei && aType>=StartRealNuclei)
            return NucleonStructure[aType-StartRealNuclei][0]+NucleonStructure[aType-StartRealNuclei][1];
        NOT_IMPLEMENTED // treat other particles when added
        return 0.;
    }

    double Particle::ElectricCharge(ParticleType aType){
        if(aType<EndLightParticle)
            return ElectricCharges[aType];
        if(aType<EndNuclei)
            return ElectricCharges[Proton]*Interactions::CNucleus::getZ(aType);
        if(aType<EndRealNuclei)
            return NucleonStructure[aType-StartRealNuclei][0];
        NOT_IMPLEMENTED // treat other particles when added
        return 0.;
    }

}
