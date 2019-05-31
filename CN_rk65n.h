/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------

  Implementation of a 6-5 Runge Kutta method with adaptive time step
  mostly taken from the book "The numerical analysis of ordinary differential
  equations - Runge-Kutta and general linear methods" by J.C. Butcher, Wiley,
  Chichester, 1987 and a free adaption to a 6 order Runge Kutta method
  of an ODE system with additive white noise

--------------------------------------------------------------------------*/  

using namespace std;

#ifndef CN_RK65N_H
#define CN_RK65N_H

#include "CN_NeuronModel.h"
#include <cmath>
#include <cfloat>

class rk65n
{
private:
  double a[9][8];
  double b[9];

  double newdt, dtx, theEps;
  double *Y[9];
  double *F[9];
  double *y5;
  double aF;
  double delta;
  int i, j, k;

protected:
  int N;
  double maxdt, eps, abseps, releps; 

public:
  rk65n(int, double, double, double, double);
  ~rk65n();
  double integrate(double *, double *, NeuronModel *, double);
  int Dim() { return N; }
};

#endif
