/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#include "CN_DCInput.h"
#include "CN_synapse.cc"

DCInput:: DCInput(neuron *target, double inI):
  synapse((neuron *) NULL, target, DCIVARNO, DCPNO, DCINPUT)
{
  p[0]= inI;               // injected current
}


DCInput::~DCInput()
{
}

void DCInput::set_I(double inI)
{
  p[0]= inI;
}

double DCInput::Isyn(double *x)
{
  return p[0];
}


