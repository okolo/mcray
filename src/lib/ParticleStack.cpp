/*
 * ParticleStack.cpp
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

#include "ParticleStack.h"
#include <signal.h>
#include "Utils.h"
#include <time.h>
#include <algorithm>

namespace mcray
{   

ParticleStack::ParticleStack(int aDebugOutputPeriod):fLastDebugOutputTime(0),fDebugOutputPeriod(aDebugOutputPeriod) {
}

ParticleStack::~ParticleStack() {
	for(std::vector<Particle*>::const_iterator  it = fPrimaryParticles.begin(); it != fPrimaryParticles.end(); it++)
		delete (*it);
	fPrimaryParticles.clear();
}


void ParticleStack::AddSecondaries(const std::vector<Particle> &aParticles)
{
#pragma omp critical (ParticleStack)
	{
	for(std::vector<Particle>::const_iterator  it = aParticles.begin(); it != aParticles.end(); it++)
		fParticles.push_back(*it);
//#ifdef _DEBUG
	time_t t = time(0);
	if(fDebugOutputPeriod>=0 && t>fLastDebugOutputTime + fDebugOutputPeriod)
	{
		std::cerr << fOutputPrefix << "particle stack: " << fParticles.size() << std::endl;
		fLastDebugOutputTime = t;
	}
//#endif
	}
}

void ParticleStack::AddPrimary(const Particle &aParticle)
{
#pragma omp critical (ParticleStack)
	{
		Particle* primary = new Particle(aParticle);
		primary->SourceParticle = 0;
		fPrimaryParticles.push_back(primary);
		fParticles.push_back(aParticle);
		//make sure SourceParticle field is valid until ParticleStack object is destroyed
		fParticles.back().SourceParticle = primary;
//#ifdef _DEBUG
		time_t t = time(0);
		if(fDebugOutputPeriod>=0 && t>fLastDebugOutputTime + fDebugOutputPeriod)
		{
			std::cerr << fOutputPrefix << "particle stack: " << fParticles.size() << std::endl;
			fLastDebugOutputTime = t;
		}
//#endif
	}
}

bool ParticleStack::Pop(Particle& aParticle)
{
	bool result = false;
#pragma omp critical (ParticleStack)
	{
	if(fParticles.size())
	{
		aParticle = fParticles.back();//*(fParticles.end()-1);
		fParticles.pop_back();//fParticles.erase(fParticles.end()-1);
//#ifdef _DEBUG
	time_t t = time(0);
	if(fDebugOutputPeriod>=0 && t>fLastDebugOutputTime + fDebugOutputPeriod)
	{
		std::cerr << fOutputPrefix << "particle stack: " << fParticles.size() << std::endl;
		fLastDebugOutputTime = t;
	}
//#endif
		result = true;
	}
	}
	return result;
}

void Result::Add(std::vector<Particle> aParticles)
{
	if(fFilter)///! filtering outside critical section
		for(int i=aParticles.size()-1; i>=0; i--)
			if (aParticles[i].Time < fT || (fFilter && !fFilter->Pass(aParticles[i])))
				aParticles.erase(aParticles.begin()+i);

#pragma omp critical (Result_Add)
	{
		for(std::vector<Particle>::const_iterator pit = aParticles.begin(); pit != aParticles.end(); pit++) {
			for (std::vector<SmartPtr<IOutput> >::iterator it = fOutputs.begin(); it != fOutputs.end(); it++) {
				(*it)->AddParticle(*pit);
			}
		}
        fParticleCounter++;//we only count number of Result::Add calls which should correspond to number of primary particles successfully treated
		if(fAutoFlushInterval>0 && fParticleCounter>=fAutoFlushInterval){//flush only when all particles from current cascade are added
			//to insure output is always selfconsistent
			fParticleCounter = 0;
			Flush(false);
		}
	}
}

void Result::Add(const Particle& aParticle)
{
	if(aParticle.Time < fT || (fFilter && !fFilter->Pass(aParticle)))///! filtering outside critical section
		return;
#pragma omp critical (Result_Add)
	{
        for(std::vector<SmartPtr<IOutput> >::iterator it = fOutputs.begin(); it != fOutputs.end(); it++)
        {
            (*it)->AddParticle(aParticle);
        }
        fParticleCounter++;
        if(fAutoFlushInterval>0 && fParticleCounter>=fAutoFlushInterval){
            fParticleCounter = 0;
            Flush(false);
        }
	}
}

void Result::SetAutoFlushInterval(uint aNumberOfPrimaryParticles){
#pragma omp critical (Result_Add)
    {
        fAutoFlushInterval = aNumberOfPrimaryParticles;
    }
}

void Result::Flush(bool aFinal)
{
#pragma omp critical (Result_Flush)
	{
		for (std::vector<SmartPtr<IOutput> >::iterator it = fOutputs.begin(); it != fOutputs.end(); it++) {
			(*it)->Flush(aFinal);
		}
	}
}

void Result::SetEndTime(const CosmoTime& aT)
{
	fT = aT;
}

cosmo_time Result::GetMinTravelTimeLeft(const Particle& aParticle) const
{
	return fT-aParticle.Time;
}

Result::~Result()
{
	Flush();
	fOutputs.clear();

	std::vector<Result *>::iterator pos = std::find(fCtrlCFlushList.begin(), fCtrlCFlushList.end(), this);
	if (pos != fCtrlCFlushList.end())
		fCtrlCFlushList.erase(pos);

}

std::vector<Result*> Result::fCtrlCFlushList;

void Result::EnableFlushOnCtrlC() {
    if (!fCtrlCFlushList.size()) {
        struct sigaction act;
        act.sa_handler = Result::SIGINT_handler;
        sigaction(SIGINT, &act, NULL);//register Ctrl-C handler
    }
    fCtrlCFlushList.push_back(this);
}

void Result::SIGINT_handler(int) {
#pragma omp critical (Result_Add)
    {
        for (std::vector<Result *>::iterator it = fCtrlCFlushList.begin(); it != fCtrlCFlushList.end(); it++) {
            (*it)->fAutoFlushInterval = 0;//disable autoflush
            (*it)->Flush(true);
            (*it)->fStopRequested = true;
        }
    }
}

}//namespace mcray

