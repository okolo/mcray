#pragma once
struct ForkSample {
  long double direction[3];
  double weight;
};
void ForkUnits(double mpc);
void *ForkField(const char *file, int mode);
void ForkValue(void *field, const long double *position, double *values);
double ForkScale(void *field);
void ForkDeleteField(void *field);
void *ForkCone(double angle, double theta, double phi, unsigned long seed);
ForkSample ForkNext(void *cone);
void ForkDeleteCone(void *cone);

double ForkVmfWeight(double theta,double phi,const long double* direction,double kappa);
