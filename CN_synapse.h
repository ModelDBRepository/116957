/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_SYNAPSE_H
#define CN_SYNAPSE_H

class synpase;

#include "CN_neuron.h"

class synapse
{
 protected:
  int idx;
  int enabled;
  
 public:
  int iVarNo;
  int pNo;
  int type;
  double *p;
  neuron *source;
  neuron *target;
  int trgiVar;
  
  synapse(neuron *, neuron *, int, int, int);
  virtual ~synapse();
  virtual void init(double *, double *);
  virtual void set_p(double *);
  virtual double gsyn()= 0;
  virtual void set_gsyn(double)= 0;
  virtual double Isyn(double *)= 0;
  void setIdx(int);
  int getIdx();
  virtual void derivative(double *, double *)= 0;
  virtual void noise(double *, double *) { }
};

#endif

