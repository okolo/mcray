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
        ParticleTypeEOF = EndNuclei //must be the last
    };

    typedef cosmo_time coord_type;

#define ParticleCustomDataSize 16
    class Particle
    {
    private:
    	inline void Reset() {
            memset(this,0,sizeof(Particle));Weight=1.; LastB = -1.; Pdir[2]=1.; PdirStart[2]=1.;
        }
        void GenerateId();
    public:
    	inline double Mass() const {
            return Mass(Type);
        }
    	inline int ElectricCharge() const {
            return ElectricCharge(Type);
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
        void PropagateFreely(cosmo_time dt);
        ParticleType Type;
        double Weight;
        CosmoTime Time;
        double Energy;
        long Ninteractions;//number of interactions in the interaction chain
        coord_type X[3];//comoving coordinates, source location: (0,0,0)
        coord_type Pdir[3];//momentum direction, initial value: (0,0,1)
        coord_type PdirStart[3];//momentum direction, initial value: (0,0,1)
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
        static const char* ParticleNames[ParticleTypeEOF];
        static int fNextCustomDataIndex;
        static unsigned long long fMaxId;
    };
}
#endif	/* PARTICLE_H */
