/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_DCINPUT_H
#define CN_DCINPUT_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define DCIVARNO 0
#define DCPNO 1

class DCInput: public synapse
{
 public:
  DCInput(neuron *, double);
  virtual ~DCInput();
  virtual void set_I(double);
  virtual double gsyn() {return 1.0;}
  virtual void set_gsyn(double) { }
  virtual double Isyn(double *);
  virtual void derivative(double *, double *) { }
};

#endif

