// SPDX-License-Identifier: MIT
#ifndef MCRAY_JET_H
#define MCRAY_JET_H
#include "Particle.h"
namespace mcray {
// Historical full-sphere vMF density multiplier, including its small-kappa
// unit-weight branch. The proposal is a uniform cone, not vMF sampling.
// Angles are radians and both direction arguments must be unit vectors.
double vMF_weight(double setup_theta,double setup_phi,const coord_type* pdir,double kappa);
// The current CLI supplies a normalized axis instead of two Euler angles.
double vMF_weight(const coord_type* axis,const coord_type* pdir,double kappa);
}
#endif
