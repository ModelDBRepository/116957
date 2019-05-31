/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_NEURONMODEL
#define CN_NEURONMODEL

#include "CN_NeuronModel.h"
#include "CN_neuron.h"
#include "CN_synapse.h"
#include "CN_neuron.cc"
#include "CN_synapse.cc"

// for inheritance ...
NeuronModel::NeuronModel()
{
}

NeuronModel::NeuronModel(list<neuron *> *ineurs, list<synapse *> *isyns,
			 int &N, ostream &msgos)
{
  int iVarCnt= 1;
  neuron *n;
  synapse *s;
  
  neurs= ineurs;
  syns= isyns;

  forall(*neurs, niter) {
    n= *niter;
    n->setIdx(iVarCnt);
    iVarCnt+= n->iVarNo;
    //    cerr << iVarCnt << endl;
  }
  forall(*syns, siter) {
    s= *siter;
    s->setIdx(iVarCnt);
    iVarCnt+= s->iVarNo;
    //    cerr << iVarCnt << endl;
  }
  N= iVarCnt;
  msgos << "# we have " << N << " variables ..." << endl;
}


NeuronModel::~NeuronModel()
{ }

void NeuronModel::derivative(double *x, double *dx)
{
  dx[0]= 1.0;
  forall(*neurs, niter) {
    (*niter)->derivative(x, dx);
  }
  forall((*syns), siter) {
    (*siter)->derivative(x, dx);
  }
}

void NeuronModel::noise(double *x, double *dx)
{
  dx[0]= 0.0;
  forall(*neurs, niter) {
    (*niter)->noise(x, dx);
  }
  forall((*syns), siter) {
    (*siter)->noise(x, dx);
  }
}

#endif
