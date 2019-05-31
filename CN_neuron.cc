//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2005-08-17
//
//--------------------------------------------------------------------------


#ifndef CN_NEURON_CC
#define CN_NEURON_CC

#include "CN_neuron.h"

neuron::neuron(int inlabel, int iniVarNo, int intype, double *inp, int inpno)
{
  label= inlabel;
  iVarNo= iniVarNo;
  type= intype;
  pno= inpno;
  if (pno > 0) {
    p= new double[pno];
    set_p(inp);
  }
  start_spiking= 0;
  spiking= 0;
  spike_time= -1.0;

  // we don't know our index number yet
  idx= 0;
  enabled= 0;
}

neuron::neuron(int inlabel, int iniVarNo, int intype, vector<int> inpos,
	       double *inp, int inpno)
{
  label= inlabel;
  iVarNo= iniVarNo;
  type= intype;
  pno= inpno;
  if (pno > 0) {
    p= new double[pno];
    set_p(inp);
  }
  pos= inpos;
  start_spiking= 0;
  spiking= 0;
  spike_time= -1.0;

 // we don't know our index number yet
  idx= 0;
  enabled= 0;
}

neuron::~neuron()
{
  forall(den, den_it) {
    (*den_it)->target= NULL;
  }
  delete[] p;
}

void neuron::set_p(double *inp)
{
  for (int i= 0; i < pno; i++) p[i]= inp[i];
}

void neuron::spike_detect(double *x)
{
  assert(enabled);
  if (E(x) >= SPK_V_THRESH)
  {
    if (!spiking)
    {
      start_spiking= 1;
      spiking= 1;
      spike_time= x[0];
    }
    else start_spiking= 0;
  }
  else {
    spiking= 0;
    start_spiking= 0;
  }
}

void neuron::init(double *x, double *iniVars)
{
  assert(enabled);
  for (int i= 0; i < iVarNo; i++)
  {
    x[idx+i]= iniVars[i];
  }
  start_spiking= 0;
  spiking= 0;
  spike_time= -1.0;
}

void neuron::setIdx(int inidx)
{
  assert(!enabled);
  idx= inidx;
  enabled= 1;
}

void neuron::noise(double *x, double *dx)
{
  for (int i= 0; i < iVarNo; i++) {
    dx[i]= 0.0;
  }
}

#endif

