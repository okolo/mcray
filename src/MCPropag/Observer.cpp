/*
 * Jet.cpp
 *
 *  Created on: Sept 15, 2021
 *      Author: dka
 */


#include <math.h>
#include <iostream>

#include "Observer.h"
#include "Particle.h"

namespace MCPropag {

using namespace std;

/*Observer::Observer(){
	this->type=-1;
	cerr<<"Observer was not inited, ignores new code"<<endl;
}//*/

Observer::Observer(int type, double x,double y,double z, double R){
	this->type=type;
	this->x=x;
	this->y=y;
	this->z=z;
	this->R=R;
	this->R2=R*R;
	this->R3=R/1.7320508075688772;
	this->Rdt=0.001*R;
	if(type==2)this->Rdt = 1.e-7*R;
	if(type==-1)cerr<<"Observer: Observer was inited as -1 type, ignores new code";//<<endl;
	else cerr<<"Observer: Observer is non trivial: x="<<this->x<<" y="<<this->y<<" z="<<this->z<<", R="<<this->R<<" (inner units)";//<<endl;
	if(type==2)cerr<<", inversed condition";
	cerr<<endl;
}

bool Observer::weInTarget(const Particle& aParticle) const{
	if(type==-1)return false; // Pseudo observer, no detection
	double dx,dy,dz;
	if(type==2){
		dz=aParticle.X[2]-z;// TODO: optimize it
		dy=aParticle.X[1]-y;
		dx=aParticle.X[0]-x;
		bool res=dx*dx+dy*dy+dz*dz>R2;
		return res;
	}
//	cerr<<"weInTarget: Particle z="<<aParticle.X[2]<<endl;
	if(fabs(dz=aParticle.X[2]-z)>R)return false;
	if(fabs(dy=aParticle.X[1]-y)>R)return false;
	if(fabs(dx=aParticle.X[0]-x)>R)return false;
	bool res=dx*dx+dy*dy+dz*dz<R2;
//	if(res)cerr<<"WE IN TARGET! ("<<aParticle.X[0]<<","<<aParticle.X[1]<<","<<aParticle.X[2]<<")"<<endl;
	return res;
}

cosmo_time Observer::patchDtIfNeed(const Particle& aParticle,cosmo_time maxDt){
	if(type==-1)return maxDt;
	double dx,dy,dz;
	dz=fabs(aParticle.X[2]-z);
	dy=fabs(aParticle.X[1]-y);
	dx=fabs(aParticle.X[0]-x);
//	cerr<<"patchDtIfNeed: (dx,dy,dz)=("<<dx<<","<<dy<<","<<dz<<")"<<endl;
//	dx=sqrt(dx*dx+dy*dy+dz*dz)-R; // TODO: Optimize it
	if(type==2){
		dx = R - sqrt(dx*dx + dy*dy + dz*dz); // TODO: Optimize it
	}else{
		if(dy>dx)dx=dy;
		if(dz>dx)dx=dz;
		dx-=R;
	}
	if(dx<=0)return Rdt;
	if(maxDt>dx)maxDt = dx;
	if(maxDt<Rdt)maxDt = Rdt; // Do not minimize to infinity
	return maxDt;
}

} /* namespace MCPropag */
