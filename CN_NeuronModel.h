/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------

The wrapper class which takes lists of pointers to neurons and synapses
which are networked to a neural system and assembles a common state
vector and handles the derivatives. At the same time it serves the neurons
and synapses their state at any given time and allows them to adjust their
parameters.

--------------------------------------------------------------------------*/  

using namespace std;

#ifndef CN_NEURONMODEL_H
#define CN_NEURONMODEL_H

#include <iostream>
#include <list>

class neuron;
class synapse;

class NeuronModel
{
 private:
  list<neuron *> *neurs;
  list<neuron *>::iterator niter;
  list<synapse *> *syns;
  list<synapse *>::iterator siter;
  
public:
  NeuronModel();
  NeuronModel(list<neuron *> *, list<synapse *> *,
	      int &, ostream &); 
  virtual ~NeuronModel();
  virtual void derivative(double *, double *);
  virtual void noise(double *, double *);
};

#endif
