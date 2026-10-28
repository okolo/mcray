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
        Nucleus = EndLightParticle,
        ParticleTypeEOF = Nucleus + 1 //must be the last
    };

    typedef cosmo_time coord_type;

#define ParticleCustomDataSize 16
    class Particle
    {
    private:
        unsigned int fA = 0, fZ = 0;
        double fNuclearMass = 0.;
    public:
        enum class NuclearMassModel { Unset, PSB, TALYS };
    private:
        NuclearMassModel fMassModel = NuclearMassModel::Unset;
        void GenerateId();
    public:
        double Mass() const;
        int ElectricCharge() const;
        unsigned int MassNumber() const;
        unsigned int AtomicNumber() const;
        void ValidateIdentity() const;
        void SetType(ParticleType type);
        void SetPsbNucleus(unsigned int a);
        void SetNucleus(unsigned int a, unsigned int z, double mass, NuclearMassModel model);
        NuclearMassModel MassModel() const { return fMassModel; }
        std::string Name() const;
        // Copies retain identity, IDs, ancestry and non-owning custom pointers.
        Particle(const Particle&) = default;
        Particle& operator=(const Particle&) = default;
        inline Particle(ParticleType aType, cosmo_time aZsource) : Type(aType) {
            Pdir[2] = 1.;
            fProductionTime.setZ(aZsource);
            Time.setZ(aZsource);
            SourceParticle = this;
            GenerateId();
        }
        Particle(unsigned int a, cosmo_time z) : Particle(Nucleus, z) { SetPsbNucleus(a); }

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
        double Weight = 1.;
        CosmoTime Time;
        double Energy = 0.;
        long Ninteractions = 0;//number of interactions in the interaction chain
        coord_type X[3]{};//comoving coordinates, source location: (0,0,0)
        coord_type Pdir[3]{};//momentum direction, initial value: (0,0,1)
        unsigned long long id = 0;
        unsigned long long fPrevId = 0;

    // TODO: interaction specific attributes should be stored separately in custom structures which can
    // be accessed via pointers stored in interactionData field
        void* interactionData[ParticleCustomDataSize]{};//used to store interaction specific attributes as pointers
        static int ReserveInteractionDataSlot();//returns index to access and store interaction data or -1 if no space left

        const Particle* SourceParticle = nullptr;
        double Deflection2 = 0.;//uncorrelated deflection angle squared
        double CorrelatedDeflection = 0.;//current value of correlated deflection
        mutable double CorrelatedBpath = 0.;// travel distance within B-field coherence length
        CosmoTime fProductionTime; // time when the particle was produced either in source or on the way as secondary from interaction;
        CosmoTime fCascadeProductionTime; // time when primary EM cascade particle was produced by hadron
        double dt = 0.; // particle time delay;
        mutable double LastB = -1.; //last transverce magnetic field within B-field coherence length

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

