/*
 * Cosmology.h
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

#ifndef COSMOLOGY_H_
#define COSMOLOGY_H_

#include "Utils.h"
#include "TableFunction.h"
#include <math.h>

namespace mcray {

using namespace Utils;
typedef long double cosmo_time;

class Cosmology
{
public:
	Cosmology();
	void Init(cosmo_time aMaxZ=100., cosmo_time aOmegaVac=0.692, cosmo_time aHubbleKm_sec_Mpc=67.8, int aAccuracy=100000);
// z -> t, dz -> dt and back transformation
	inline cosmo_time z2t0(cosmo_time _z) const;//t0==0 at z=0
	inline cosmo_time t02z(cosmo_time _t) const;//t0==0 at z=0
	inline cosmo_time z2t(cosmo_time _z) const;//t==0 at z=infinity
	inline cosmo_time t2z(cosmo_time _t) const;//t==0 at z=infinity
	inline cosmo_time z2d(cosmo_time _z) const;
	inline cosmo_time d2z(cosmo_time _d) const;
	inline cosmo_time getLv() const{return Lv;};
	inline cosmo_time getAgeOfUniverse() const{return t0;};
	void SuppressionTest() const;
	//returns energy loss rate -1/E dE/dt due to redshift at
	inline cosmo_time eLossRate(cosmo_time aZ) const;

	// -dt/dz
	inline cosmo_time diffT(cosmo_time aZ) const;

	//Hubble parameter at given redshift
	inline cosmo_time Hubble(cosmo_time aZ) const;
	//test output
	void print() const;
	static void UnitTest();
	inline bool IsInitialized() { return t0!=0; }
	inline cosmo_time RelError() const { return 1./Accuracy; }
private:
	//test
	static cosmo_time Hy2F1(cosmo_time a, cosmo_time b, cosmo_time c, cosmo_time z);
	cosmo_time TestIntegralF(cosmo_time r) const;
	cosmo_time TestIntegral(cosmo_time r) const;

	inline cosmo_time diffD(cosmo_time a) const;
	void initZ2D();
	void initZ2T();

	cosmo_time H;//Hubble parameter at z=0 (internal units)
	cosmo_time Lv;// omega_lambda at z=0
	cosmo_time Lm;// omega_matter at z=0
	cosmo_time sqrtLv;//sqrt(Lv)
	cosmo_time sqrtLm;//sqrt(Lv)
	cosmo_time t0;//universe age at z=0
	cosmo_time fMaxZ;
	int Accuracy;
	std::vector<cosmo_time>	zt;
	std::vector<cosmo_time>	tt;
	std::vector<cosmo_time>	z;
	std::vector<cosmo_time>	d;
	SafePtr< LinearFuncX<cosmo_time> >	z2dFunc;
	SafePtr< LinearFuncX<cosmo_time> >	d2zFunc;
	SafePtr< LinearFuncX<cosmo_time> >	z2tFunc;
	SafePtr< LinearFuncX<cosmo_time> >	t2zFunc;
};

inline cosmo_time Cosmology::t2z(cosmo_time _t) const
{
	 return t02z(_t-t0);
}

inline cosmo_time Cosmology::z2t0(cosmo_time _z) const//t0==0 at z=0
{
	ASSERT(z2tFunc->InTableRange(_z));
	return -z2tFunc->f(_z);
}

inline cosmo_time Cosmology::t02z(cosmo_time _t0) const//t0==0 at z=0
{
	cosmo_time timeAgoFromNow = -_t0;
	if(timeAgoFromNow<0)
	{
		ASSERT(-timeAgoFromNow/t0<1e-12);//roundoff error
		timeAgoFromNow=0.;
	}
	ASSERT(t2zFunc->InTableRange(timeAgoFromNow));
	return t2zFunc->f(timeAgoFromNow);
}

inline cosmo_time Cosmology::z2t(cosmo_time _z) const
{
	ASSERT(z2tFunc->InTableRange(_z));
	return t0 - z2tFunc->f(_z);
}

inline cosmo_time Cosmology::eLossRate(cosmo_time aZ) const
{
	return 1./(1.+aZ)/diffT(aZ);
}

inline cosmo_time Cosmology::diffD(cosmo_time a) const
{
	return 1./( H*sqrt(a*(Lm+Lv*a*a*a)) );
}

inline cosmo_time Cosmology::z2d(cosmo_time _z) const
{
	ASSERT(z2dFunc->InTableRange(_z));
	return z2dFunc->f(_z);
}

inline cosmo_time Cosmology::d2z(cosmo_time _d) const
{
	ASSERT(d2zFunc->InTableRange(_d));
	return d2zFunc->f(_d);
}

inline cosmo_time Cosmology::diffT(cosmo_time aZ) const
{
	cosmo_time z1 = 1.+aZ;
	return 1./H/z1/sqrt(Lm*z1*z1*z1+Lv);
}

inline cosmo_time Cosmology::Hubble(cosmo_time aZ) const
{
	return 1./((1.+aZ)*diffT(aZ));
}

extern Cosmology cosmology;

class CosmoTime
{
public:
	//construct with z=0
	CosmoTime():m_z(0.),m_t(0.){}

	//constructor with redshift
	CosmoTime(cosmo_time aRedshift) {setZ(aRedshift);}

	//copy constructor
	CosmoTime(const CosmoTime& aTime):m_z(aTime.m_z), m_t(aTime.m_t){}

	inline operator cosmo_time() const { return t(); }
	inline CosmoTime& operator=(const CosmoTime& aTime) { m_z = aTime.m_z; m_t = aTime.m_t; return *this; }
	inline CosmoTime& operator=(cosmo_time aT) { setT(aT); return *this; }
	inline CosmoTime& operator += (cosmo_time aDt)
	{
		cosmo_time curT = t();
		cosmo_time newT = curT+aDt;
		if(newT>0){
			ASSERT(newT/aDt<cosmology.RelError());
			newT = 0.;//fix roundoff error
		}
		setT(newT); return *this;
	}

	inline bool operator<(const CosmoTime& aT) const
	{
		ASSERT(m_z>=0. || m_t<=0);
		if(m_z<0 || aT.m_z<0)
			return t()<aT.t();
		else
			return aT.m_z<m_z;
	}

	//redshift value
	inline cosmo_time z() const
	{
		ASSERT(m_z>=0. || m_t<=0);
		if(m_z<0)
			m_z = cosmology.t02z(m_t);//lazy initialization
		return m_z;
	}

	//time difference from z=0 (<0 for past moments)
	inline cosmo_time t() const
	{
		ASSERT(m_z>=0. || m_t<=0);
		if(m_t>0.)
			m_t = cosmology.z2t0(m_z);//lazy initialization
		return m_t;
	}

	//set redshift value z
	inline void setZ(cosmo_time aZ)
	{
		ASSERT(aZ>=0.);
		m_z = aZ;
		m_t = 1;//enable lazy initialization
	}

	//set redshift using time difference from z=0
	inline void setT(cosmo_time aT)
	{
		ASSERT(aT<=0.);
		m_z = -1;//enable lazy initialization
		m_t = aT;
	}

	std::string ToString() const {
		std::ostringstream logStr;
		logStr << "z=" << m_z << ", t=" <<m_t;
		return logStr.str();
	}
private:
	mutable cosmo_time m_z;
	mutable cosmo_time m_t;//time difference between z=0 and this z (always < 0)
};

} /* namespace mcray */
#endif /* COSMOLOGY_H_ */
