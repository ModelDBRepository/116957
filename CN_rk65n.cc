/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

using namespace std;

#ifndef CN_RK65N_CC
#define CN_RK65N_CC

#include "CN_rk65n.h"
#include <iostream>

#define min(X,Y) (X<Y? X:Y)
#define max(X,Y) (X>Y? X:Y)

#define sixth 0.166666666666667



rk65n::rk65n(int dim, double inmaxdt, double ineps, double inabseps,
	     double inreleps)
{
  N= dim;
  maxdt= inmaxdt;
  eps= ineps;
  abseps= inabseps;
  releps= inreleps;
  for (int i= 0; i < 9; i++) {
    Y[i]= new double[N];
    F[i]= new double[N];
    for (int j= i; j < 8; j++) {
      a[i][j]= 0.0;
    }
  }
  y5= new double[N];
  
  a[1][0]= 0.111111111111111;
  a[2][0]= 0.555555555555556e-1;
  a[2][1]= 0.555555555555556e-1;
  a[3][0]= 0.416666666666667e-1;
  a[3][1]= 0.0;
  a[3][2]= 0.125;
  a[4][0]= 0.166666666666667;
  a[4][1]= 0.0;
  a[4][2]= -0.5;
  a[4][3]= 0.666666666666667;
  a[5][0]= 0.1875e+1;
  a[5][1]= 0.0;
  a[5][2]= -0.7875e+1;
  a[5][3]= 0.7e+1;
  a[5][4]= -0.5;
  a[6][0]= -0.4227272727272727e+1;
  a[6][1]= 0.0;
  a[6][2]= 0.176995738636364e+2;
  a[6][3]= -0.142883522727273e+2;
  a[6][4]= 0.522017045454545;
  a[6][5]= 0.104403409090909e+1;
  a[7][0]= 0.840622673179752e+1;
  a[7][1]= 0.0;
  a[7][2]= -0.337303717185049e+2;
  a[7][3]= 0.271460231129622e+2;
  a[7][4]= 0.342046929709216;
  a[7][5]= -0.184653767923258e+1;
  a[7][6]= 0.577349465373733;
  a[8][0]= 0.128104575163399;
  a[8][1]= 0.0;
  a[8][2]= 0.0;
  a[8][3]= -0.108433734939759;
  a[8][4]= 0.669375;
  a[8][5]= -0.146666666666667;
  a[8][6]= 0.284444444444444;
  a[8][7]= 0.173176381998583;

  b[0]= 0.567119155354449e-1;
  b[1]= 0.0;
  b[2]= 0.0;
  b[3]= 0.210909572355356;
  b[4]= 0.141490384615385;
  b[5]= 0.202051282051282;
  b[6]= 0.253186813186813;
  b[7]= 0.843679809736684e-1;
  b[8]= 0.512820512820513e-1;
}

rk65n::~rk65n()
{
  for (int i= 0; i < 9; i++) {
    delete[] Y[i];
    delete[] F[i];
  }
  delete[] y5;
}

double rk65n::integrate(double *y, double *yn,
			NeuronModel *model, double dt)
{
  // calculate iterative terms rk65_Y[__i] and rk65_F[__i] (to sixth order)
  for (i= 0; i < 9; i++)
  {
    for (k= 0; k < N; k++)
    {
      aF= 0.0;
      for (j= 0; j < i; j++)
	aF+= a[i][j]*F[j][k];
      Y[i][k]= y[k]+dt*aF;
    }
    model->derivative(Y[i], F[i]);
  }

  // sum up rk65_Y[__i] and rk65_F[__i] to build 5th order scheme -> rk65_y5
  for (k= 0; k < N; k++)
  {
    aF= 0.0;
    for (j= 0; j < 8; j++) aF+= a[8][j]*F[j][k];
    y5[k]= y[k]+ dt*aF;
  }

  // sum up rk65_Y[__i] and rk65_F[__i] to build 6th order scheme -> yn
  for (k= 0; k < N; k++)
  {
    aF= 0.0;
    for (j= 0; j < 9; j++) aF+= b[j]*F[j][k];
    yn[k]= y[k]+ dt*aF;
  }

  // determine minimal necessary new dt to get error < theEps based on the
  // difference between results rk65_y5 and yn
  dtx= maxdt;
#ifdef DEBUG
  int min_var= -1;
#endif
  for (k= 0; k < N; k++)
  {
    theEps= max(abseps, min(eps, fabs(releps*yn[k])));
    delta= abs(yn[k]-y5[k]);
    if (delta > DBL_MIN) {
      newdt= exp(sixth*(log(theEps)-log(delta)))*dt;
      if (newdt < dtx) {
#ifdef DEBUG
	min_var= k;
#endif
	dtx= newdt;
      }
    }
  }
#ifdef DEBUG
  cerr << min_var << " " << dtx << endl;
#endif
  return dtx;
}

#undef min
#undef max

#endif
