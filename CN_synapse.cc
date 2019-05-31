/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_SYNAPSE_CC
#define CN_SYNAPSE_CC

#include "CN_synapse.h"

synapse::synapse(neuron *insource, neuron *intarget, int iniVarNo,
		 int inpNo, int intype)
{
  assert(intarget != NULL);
  source= insource;
  target= intarget;
  target->den.push_back(this);
  iVarNo= iniVarNo;
  pNo= inpNo;
  if (pNo > 0) p= new double[pNo];
  type= intype;

  // we don't know our index number yet
  idx= 0;
  enabled= 0;
}

synapse::~synapse()
{
  if (pNo > 0) delete[] p;
}

void synapse::set_p(double *inp)
{
  for (int i= 0; i < pNo; i++) p[i]= inp[i];
}

void synapse::init(double *x, double *iniVars)
{
  assert(enabled);
  for (int i= 0; i < iVarNo; i++) {
    x[idx+i]= iniVars[i];
  }
}


void synapse::setIdx(int inidx)
{
  assert(!enabled);
  idx= inidx;
  enabled= 1;
}

int synapse::getIdx()
{
  return idx;
}

#endif
