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


//--------------------------------------------------------------------------
// General abstract class for a neuron
//--------------------------------------------------------------------------

#ifndef CN_NEURON_H
#define CN_NEURON_H

double SPK_V_THRESH= 0.0;

class neuron;
class synapse;

#include "CN_base.h"
#include "CN_synapse.h"

class neuron
{
 protected:
  int enabled;
  
 public:
  int idx;
  int label;
  int iVarNo;
  int type;
  int pno;
  double *p;
  vector<int> pos;
  int spiking, start_spiking;
  double spike_time;

  list<synapse *> den;       // list of all incoming synapses
  list<synapse *>::iterator den_it;
  
  neuron(int, int, int, double *, int);
  neuron(int, int, int, vector<int>, double *, int);
  virtual ~neuron();
  virtual void set_p(double *);
  virtual void init(double *, double *);
  virtual void spike_detect(double *);
  void setIdx(int);
  int getIdx() { return idx;}
  virtual double E(double *)= 0;
  virtual void derivative(double *, double *)= 0;
  virtual void noise(double *, double *);
};

#endif
