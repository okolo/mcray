/*
 * Particle.h
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

#ifndef PARTICLE_H
#define	PARTICLE_H

#include "Cosmology.h"
#include <cstring>
#include "Units.h"
#include <math.h>

namespace mcray
{
    enum ParticleType
    {// order of masses in the array Particle::massesMeV must follow the order of particles below
        Electron = 0,
        Positron,
        Photon,
    	NeutrinoE,
    	NeutrinoM,
    	NeutrinoT,
    	NeutrinoAE,
    	NeutrinoAM,
    	NeutrinoAT,
    	Neutron,
    	Proton,
        EndLightParticle,
        StartNuclei = EndLightParticle,
        A02 = StartNuclei,
        A03,
        A04,
        A05,
        A06,
        A07,
        A08,
        A09,
        A10,
        A11,
        A12,
        A13,
        A14,
        A15,
        A16,
        A17,
        A18,
        A19,
        A20,
        A21,
        A22,
        A23,
        A24,
        A25,
        A26,
        A27,
        A28,
        A29,
        A30,
        A31,
        A32,
        A33,
        A34,
        A35,
        A36,
        A37,
        A38,
        A39,
        A40,
        A41,
        A42,
        A43,
        A44,
        A45,
        A46,
        A47,
        A48,
        A49,
        A50,
        A51,
        A52,
        A53,
        A54,
        A55,
        A56,
        EndNuclei,
        StartRealNuclei=EndNuclei,
        H2=StartRealNuclei,
        H3,
        He3,
        He4,
        Li6,
        Li7,
        Be7,
        Be9,
        Be10,
        Be11,
        B10,
        B11,
        C10,
        C11,
        C12,
        C13,
        C14,
        C15,
        C16,
        N13,
        N14,
        N15,
        N16,
        N17,
        O14,
        O15,
        O16,
        O17,
        O18,
        O19,
        O20,
        O21,
        O22,
        F17,
        F18,
        F19,
        F20,
        F21,
        F22,
        F23,
        Ne18,
        Ne19,
        Ne20,
        Ne21,
        Ne22,
        Ne23,
        Ne24,
        Na21,
        Na22,
        Na23,
        Na24,
        Na25,
        Mg22,
        Mg23,
        Mg24,
        Mg25,
        Mg26,
        Mg27,
        Mg28,
        Al24,
        Al25,
        Al26,
        Al27,
        Al28,
        Al29,
        Al30,
        Si26,
        Si27,
        Si28,
        Si29,
        Si30,
        Si31,
        Si32,
        Si33,
        Si34,
        P29,
        P30,
        P31,
        P32,
        P33,
        P34,
        P35,
        P36,
        P37,
        S31,
        S32,
        S33,
        S34,
        S35,
        S36,
        S37,
        S38,
        S39,
        S40,
        S41,
        Cl33,
        Cl34,
        Cl35,
        Cl36,
        Cl37,
        Cl38,
        Cl39,
        Cl40,
        Cl41,
        Cl42,
        Cl43,
        Ar35,
        Ar36,
        Ar37,
        Ar38,
        Ar39,
        Ar40,
        Ar41,
        Ar42,
        Ar43,
        Ar44,
        Ar45,
        Ar46,
        K38,
        K39,
        K40,
        K41,
        K42,
        K43,
        K44,
        K45,
        K46,
        K47,
        K48,
        K49,
        Ca40,
        Ca41,
        Ca42,
        Ca43,
        Ca44,
        Ca45,
        Ca46,
        Ca47,
        Ca48,
        Ca49,
        Ca50,
        Sc43,
        Sc44,
        Sc45,
        Sc46,
        Sc47,
        Sc48,
        Sc49,
        Sc50,
        Sc51,
        Ti44,
        Ti45,
        Ti46,
        Ti47,
        Ti48,
        Ti49,
        Ti50,
        Ti51,
        Ti52,
        V47,
        V48,
        V49,
        V50,
        V51,
        V52,
        V53,
        Cr48,
        Cr49,
        Cr50,
        Cr51,
        Cr52,
        Cr53,
        Cr54,
        Mn51,
        Mn52,
        Mn53,
        Mn54,
        Mn55,
        Fe52,
        Fe53,
        Fe54,
        Fe55,
        Fe56,
        Fe57,
        Fe58,
        Fe59,
        Fe60,
        Fe61,   // 253
        EndRealNuclei, // 254
        ParticleTypeEOF = EndRealNuclei //must be the last
    };

    typedef cosmo_time coord_type;

#define ParticleCustomDataSize 16
    class Particle
    {
    private:
    	inline void Reset() {
            memset(this,0,sizeof(Particle));Weight=1.; LastB = -1.; Pdir[2]=1.;
        }
        void GenerateId();
    public:
    	inline double Mass() const {
            return Mass(Type);
        }
    	inline int ElectricCharge() const {
            return ElectricCharge(Type);
        }
        inline int AtomicMass() const {
            return AtomicMass(Type);
        }
        //Construct primary particle
    	inline Particle(ParticleType aType, cosmo_time aZsource){
            Reset(); Type=aType; fProductionTime.setZ(aZsource) ; Time.setZ(aZsource); SourceParticle = this;
            GenerateId();
        }

    	inline Particle& operator=(const Particle& aParticle) {
            memcpy(this, &aParticle, sizeof(Particle)); return *this;
        }
    	inline coord_type beta() const {
            coord_type m=Mass();
            coord_type gamma_1 = m/Energy;
            coord_type beta2 = 1. - gamma_1*gamma_1;
            return sqrt(beta2);
        }

        inline void SetParent(const Particle& aParticle){
            GenerateId();
            if(id != aParticle.id){
                fPrevId = aParticle.id;
                fProductionTime = Time;
            }
        }

    	static double Mass(ParticleType aType);
    	static double ElectricCharge(ParticleType aType);
        static double AtomicMass(ParticleType aType);
        void PropagateFreely(cosmo_time dt);
        ParticleType Type;
        double Weight;
        CosmoTime Time;
        double Energy;
        long Ninteractions;//number of interactions in the interaction chain
        coord_type X[3];//comoving coordinates, source location: (0,0,0)
        coord_type Pdir[3];//momentum direction, initial value: (0,0,1)
        unsigned long long id;
        unsigned long long fPrevId;

    // TODO: interaction specific attributes should be stored separately in custom structures which can
    // be accessed via pointers stored in interactionData field
        void* interactionData[ParticleCustomDataSize];//used to store interaction specific attributes as pointers
        static int ReserveInteractionDataSlot();//returns index to access and store interaction data or -1 if no space left

        const Particle* SourceParticle;
        double Deflection2;//uncorrelated deflection angle squared
        double CorrelatedDeflection;//current value of correlated deflection
        mutable double CorrelatedBpath;// travel distance within B-field coherence length
        CosmoTime fProductionTime; // time when the particle was produced either in source or on the way as secondary from interaction;
        CosmoTime fCascadeProductionTime; // time when primary EM cascade particle was produced by hadron
        double dt; // particle time delay;
        mutable double LastB; //last transverce magnetic field within B-field coherence length

        inline double JetOpenningAngle() const
        {
        	return DeflectionAngle()-ObservationAngle();
        }

        inline CosmoTime LastDeflectionTime() const
        {
        	return ElectricCharge() ? Time : fProductionTime;
        }

        inline double ObservationAngle() const
        {
        	double beta = DeflectionAngle();
        	if(beta==0.)
        		return 0.;
            double sinBeta=sin(beta);
            double pathToLastDeflection = LastDeflectionTime().t()-SourceParticle->Time.t();//distance from the source for the parent electron/positron (valid for small z source and deflection)
            double distanceToSource = Time.t()-SourceParticle->Time.t();// distance between the source and the observer  (valid for small z source)
            return asin(pathToLastDeflection/distanceToSource*sinBeta);// observation angle
        }

        inline double DeflectionAngle() const
        {
        	return sqrt(Deflection2 + CorrelatedDeflection*CorrelatedDeflection);
        }

        inline double TimeDelay() const
        {
        	if(ElectricCharge())
        		return 0.;//the approximation used below works only for neutral particles ( LastDeflectionTime() > Time.t() )
            double sinBeta=sin(DeflectionAngle());
            double pathToLastDeflection = LastDeflectionTime().t()-SourceParticle->Time.t();//xx - distance from the source for the parent electron/positron (valid for small z source and deflection)
            double distanceToSource = Time.t()-SourceParticle->Time.t();// distance between the source and the observer  (valid for small z source)
            return 2*pathToLastDeflection*(1.0-pathToLastDeflection/distanceToSource)*sinBeta*sinBeta;//  time delay
        }
        inline static const char* Name(ParticleType aType)
        {
            return ParticleNames[aType];
        }
        std::string ToString();
    private:
        static const double MassesMeV[EndLightParticle];
        static const int ElectricCharges[EndLightParticle];
        static const int NucleonStructure[EndRealNuclei-StartRealNuclei][2];
        static const char* ParticleNames[ParticleTypeEOF];
        static int fNextCustomDataIndex;
        static unsigned long long fMaxId;
    };
}
#endif	/* PARTICLE_H */

