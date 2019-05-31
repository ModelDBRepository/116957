/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-01-25
  
--------------------------------------------------------------------------*/

using namespace std;

#include <sstream>
#include <fstream>
#include "CN_LPneuronNT.h"
#include "CN_DCInput.h"
#include "CN_NeuronModel.h"
#include "CN_rk65n.h"

#include "CN_LPneuronNT.cc"
#include "CN_DCInput.cc"
#include "CN_NeuronModel.cc"
#include "CN_rk65n.cc"

double rk65_MINDT= 0.05;
double rk65_eps= 1e-12;
double rk65_relEps= 1e-9;
double rk65_absEps= 1e-16;
double mindt= 1e-6;

int main(int argc, char *argv[])
{
  if (argc != 5) {
    cerr << "usage: testNT <input_current> <delta I> <step no> <basename>" << endl;
    exit(1);
  }

  cerr << "call was: ";
  for (int i= 0; i < argc; i++) {
    cerr << argv[i] << " ";
  }
  cerr << endl;
		 
  double inI= atof(argv[1]);
  double dI= atof(argv[2]);
  int steps= atoi(argv[3]);
  list<neuron *> neurs;
  list<synapse *> syns;
  
  vector<int> dummy(3);
  LPRneuron n(1, dummy);
  neurs.push_back(&n);
    
  n.set_p(stdLPR_p);
  DCInput s(&n, inI);
  double told, t, dt, dtx;
  double theI;
  stringstream name;
  char thename[80];
  ofstream os;
  os.precision(10);
  
  for (int i= 0; i < LPR_PNO; i++) {
    cerr << "# " << LPR_p_text[i] << " ";
    cerr << n.p[i] << endl;
  }

  double *x, *xn, *tmp;
  int N;
  NeuronModel model(&neurs, &syns, N, cerr);
  x= new double[N];
  xn= new double[N];
  rk65n machine(N, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);

  cout.precision(10);
  for (int k= 0; k < steps; k++) {
    theI= inI+k*dI;

    name.clear();
    name << argv[4] << "." << theI << ".dat" << ends;
    name >> thename;
    os.open(thename);
    
    s.set_I(theI-2.0);
    told= -0.1;
    t= 0.0;
    dt= 0.0001;
    dtx= 0.005;
    n.init(x, LPR_INIVARS);
    
    x[0]= 0;
    while (x[0] < 20000)
    {
      if ((x[0] > 5000) && (x[0] - told) >= 0.1)
      {
	os << x[0] << " " << n.E(x) << endl;
	told= x[0];
      }
      if (isnan(n.E(x)))  {
	cerr << "nan encountered!" << endl;
	exit(1);
      }
      
      dtx= machine.integrate(x, xn, &model, dt);
      dtx= min(dtx, 2.0*dt);

      tmp= x; x= xn; xn= tmp;
      dt= dtx;
    }
    os.close();
  }
  delete[] x;
  delete[] xn;
  
  return 0;
}


