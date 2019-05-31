/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-01-25
  
--------------------------------------------------------------------------*/

#ifndef LPRNEURON_CC
#define LPRNEURON_CC

#include "CN_neuron.cc"

#define efunc(X,Y,Z) (1.0/(1.0+exp(((X)-(Y))/(Z))))
#define Eaxon x[idx]
#define Esoma x[idx+15]

#define gNa p[0]
#define ENa p[1] 
#define gCa1 p[2]
#define gCa2 p[3]
#define goCa p[4]
#define EK p[5]
#define gd p[6]
#define gA p[7]
#define gh p[8]
#define EIh p[9]
#define gleak p[10]
#define Eleak p[11]
#define Cmem p[12]
#define c_iCa p[13]
#define k_Ca p[14]
#define kfac p[15]
#define Vshift p[16]
#define IDC p[17]
#define kfacfast p[18]
#define gVV p[19]
#define Cs p[20]
#define gleaks p[21]
#define RTF p[22]
#define AreaAxon p[23]
#define Area p[24]

// ICa from my own fit
#define k_aCa1 p[25]
#define k_bCa1 p[26]
#define k_aCa2 p[27]
#define V_aCa1 p[28]
#define V_bCa1 p[29]
#define V_aCa2 p[30]
#define s_aCa1 p[31]
#define s_bCa1 p[32]
#define s_aCa2 p[33]
#define P_Ca p[34]
#define Ca_out p[35]
#define V_kbCa1 p[36]
#define s_kbCa1 p[37]

#define k_oa p[38]
#define k_ob p[39]
#define V_ao1 p[40]
#define V_ao2 p[41]
#define s_ao1 p[42]
#define s_ao2 p[43]
#define f p[44]
#define c1 p[45]
#define c2 p[46]
#define c3 p[47]
#define Ca_0 p[48]

#define k_aA p[49]
#define k_bA1 p[50]
#define c_A2 p[51]
#define V_aA p[52]
#define V_bA p[53]
#define V_kbA2 p[54]
#define V_x p[55]
#define s_aA p[56]
#define s_bA p[57]
#define s_kbA2 p[58]
#define s_x p[59]

// Ih from my own fit
#define c_r p[60]
#define V_r p[61]
#define V_kr p[62]
#define s_r p[63]
#define s_kr p[64]

#define I_scale p[65]

#define E_M p[66]
#define g_M p[67]
#define k_M p[68]
#define V_M p[69]
#define s_M p[70]
#define V_kM p[71]
#define s_kM p[72]

LPRneuron::LPRneuron(int inlabel, double *the_p= LPR_p):
  neuron(inlabel, LPR_IVARNO, LPRNEURON, the_p, LPR_PNO)
{
}

LPRneuron::LPRneuron(int inlabel, vector<int> inpos, double *the_p= LPR_p):
  neuron(inlabel, LPR_IVARNO, LPRNEURON, inpos, the_p, LPR_PNO)
{
}

inline double LPRneuron::E(double *x)
{
  return x[idx+15]+Vshift;
}

void LPRneuron::currents(ostream &os, double *x)
{
  static double INa, ICa, IoCa, Id, IA, Ih, Il, IM, IVV, sum;

  // differential eqn for the axon membrane potential
  INa= pw3(x[idx+1])*x[idx+2]*gNa*(Eaxon-ENa)*AreaAxon;
  Id= pw4(x[idx+8])*gd*(Esoma-EK)*AreaAxon;
  IM= g_M*x[idx+14]*(Eaxon-E_M)*AreaAxon;
  Il= gleak*(Eaxon-Eleak)*AreaAxon;
  IVV= gVV*(Esoma-Eaxon);
  
  os << INa << " ";
  os << Id << " ";
  os << IM << " ";
  os << Il << " ";
  sum= INa + Id + IM + Il - IVV;
  os << sum << " ";
  
  ICa= (x[idx+3]*x[idx+4]*gCa1+x[idx+5]*gCa2)*
    P_Ca*Esoma*(x[idx+13]*exp(Esoma/RTF)-Ca_out)/(exp(Esoma/RTF)-1.0)*Area;
  IoCa= x[idx+6]*x[idx+7]*goCa*(Esoma-EK)*Area;
  IA= pw3(x[idx+9])*gA*(Esoma-EK)*
    (efunc(Esoma,V_x,s_x)*x[idx+10]+(1.0-efunc(Esoma,V_x,s_x))*x[idx+11])*Area;
  Ih= x[idx+12]*gh*(Esoma-EIh)*Area;
  Il= gleaks*(Esoma-Eleak)*Area;

  sum= ICa + IoCa + IA + Ih + IVV + Il;
  os << ICa << " ";
  os << IoCa << " ";
  os << IA << " ";
  os << Ih << " ";
  os << IVV << " ";
  os << Il << " ";
  os << sum << endl;
}

void LPRneuron::derivative(double *x, double *dx)
{
  static double a, b, hinf, kh, minf, km;
  static double INa, Id, ICa, IoCa, IA, Ih, Il, IM, IVV;
  static double Isyn;
  
  Isyn= 0.0;
  forall(den, den_it) {
    Isyn+= (*den_it)->Isyn(x);
  }
  
  // differential eqn for the axon membrane potential
  INa= pw3(x[idx+1])*x[idx+2]*gNa*(Eaxon-ENa)*AreaAxon;
  Id= pw4(x[idx+8])*gd*(Eaxon-EK)*AreaAxon;
  IM= g_M*x[idx+14]*(Eaxon-E_M)*AreaAxon;
  Il= gleak*(Eaxon-Eleak)*AreaAxon;
  IVV= gVV*(Esoma-Eaxon);
  
  dx[idx]= (-INa - Id - IM - Il + IVV)/Cmem; 

  // soma compartment
  a= exp(Esoma/RTF);
  ICa= (x[idx+3]*x[idx+4]*gCa1+x[idx+5]*gCa2)*P_Ca*Esoma*(x[idx+13]*a-Ca_out)/(a-1.0)*Area;
  
  IoCa= x[idx+6]*x[idx+7]*goCa*(Esoma-EK)*Area;

  a= efunc(Esoma,V_x,s_x);
  IA= pw3(x[idx+9])*gA*(Esoma-EK)*(a*x[idx+10]+(1.0-a)*x[idx+11])*Area;

  Ih= x[idx+12]*gh*(Esoma-EIh)*Area;

  Il= gleaks*(Esoma-Eleak)*Area;

  dx[idx+15]= (-ICa - IoCa - IA - Ih - Il - IVV + (IDC+Isyn)*I_scale)/Cs;

  // differential eqn for mNa
  a= (-16.64-0.32*x[idx]) / (exp(-13.0-x[idx]/4.0)-1.0);
  b= (0.28*x[idx]+7.0) / (exp(x[idx]/5.0+5.0)-1.0);
  dx[idx+1]= a*(1.0-x[idx+1])-b*x[idx+1];
  // differential eqn for hNa
  a= 0.128*exp(-2.66666666667-x[idx]/18.0);
  b= 4.0 / (exp(-5.0-x[idx]/5.0)+1.0);
  dx[idx+2]= a*(1.0-x[idx+2])-b*x[idx+2];

  // differential eqn for md
  a=  (-1.6-0.032*x[idx]) / (exp(-10.0-x[idx]/5.0)-1.0);
  b= 0.5*exp(-1.375-x[idx]/40.0);
  dx[idx+8]= a*(1.0-x[idx+8])-b*x[idx+8];

  // differential eqn for mCa1
  minf= efunc(Esoma,V_aCa1,s_aCa1);
  km= k_aCa1*kfac;
  dx[idx+3]= (minf-x[idx+3])*km;
  // differential eqn for hCa1
  hinf= efunc(Esoma,V_bCa1,s_bCa1);
  kh= k_bCa1*efunc(Esoma,V_kbCa1,s_kbCa1)*kfac;
  dx[idx+4]= (hinf-x[idx+4])*kh;

  // differential eqn for mCa2
  minf= efunc(Esoma,V_aCa2,s_aCa2);
  km= k_aCa2*kfac;
  dx[idx+5]= (minf-x[idx+5])*km;

  // differential eqn for moCa
  a= V_ao1-f*x[idx+13];
  b= V_ao2-f*x[idx+13];
  minf= efunc(Esoma,a,s_ao1)*efunc(Esoma,b,s_ao2)*(x[idx+13]/(c1+x[idx+13]));
  km= k_oa*kfac;
  dx[idx+6]= (minf-x[idx+6])*km;
  // differential eqn for hoCa
  hinf= c2/(c3+x[idx+13]);
  kh= k_ob*kfac;
  dx[idx+7]= (hinf-x[idx+7])*kh;

  // differential eqn for mA
  minf= efunc(Esoma,V_aA,s_aA);
  km= k_aA*kfac;
  dx[idx+9]= (minf-x[idx+9])*km;
  // differential eqn for hA1
  hinf= efunc(Esoma,V_bA,s_bA);
  kh= k_bA1*kfac;
  dx[idx+10]= (hinf-x[idx+10])*kh;
  // differential eqn for hA2
  kh= c_A2*efunc(Esoma,V_kbA2,s_kbA2)*kfac;
  dx[idx+11]= (hinf-x[idx+11])*kh;
    
  // differential eqn for mh
  minf= efunc(Esoma,V_r,s_r); 
  km= c_r/efunc(Esoma,V_kr,s_kr)*kfac;
  dx[idx+12]= (minf-x[idx+12])*km;
     
  // differential eqn for Ca concentration 
  dx[idx+13]= -c_iCa*ICa - k_Ca*(x[idx+13]-Ca_0);

  // differential equn for M current activation var
  minf= efunc(Eaxon,V_M,s_M);
  km= k_M*efunc(Eaxon,V_kM,s_kM)*kfac;
  dx[idx+14]= (minf-x[idx+14])*km;
}

#undef efunc
#undef Eaxon
#undef Esoma

#undef gNa 
#undef ENa 
#undef gCa1 
#undef gCa2 
#undef goCa 
#undef EK 
#undef gd 
#undef gA 
#undef gh 
#undef EIh 
#undef gleak 
#undef Eleak 
#undef Cmem 
#undef c_iCa 
#undef k_Ca 
#undef kfac 
#undef Vshift 
#undef IDC 
#undef kfacfast 
#undef gVV 
#undef Cs 
#undef gleaks 
#undef RTF 
#undef Area
#undef AreaAxon

// ICa from my own fit
#undef k_aCa1 
#undef k_bCa1 
#undef k_aCa2 
#undef V_aCa1 
#undef V_bCa1 
#undef V_aCa2 
#undef s_aCa1 
#undef s_bCa1 
#undef s_aCa2 
#undef P_Ca 
#undef Ca_out 
#undef V_kbCa1 
#undef s_kbCa1 

#undef k_oa 
#undef k_ob 
#undef V_ao1 
#undef V_ao2 
#undef s_ao1 
#undef s_ao2 
#undef f 
#undef c1 
#undef c2 
#undef c3 
#undef Ca_0 

#undef c_n 
#undef V_n 
#undef V_kn 
#undef s_n 
#undef s_kn 

#undef k_aA 
#undef k_bA1 
#undef c_A2 
#undef V_aA 
#undef V_bA 
#undef V_kbA2 
#undef V_x 
#undef s_aA 
#undef s_bA 
#undef s_kbA2 
#undef s_x 

// Ih from my own fit
#undef c_r 
#undef V_r 
#undef V_kr 
#undef s_r 
#undef s_kr 

#undef I_scale 

#undef E_M 
#undef g_M 
#undef k_M 
#undef V_M 
#undef s_M 
#undef V_kM
#undef s_kM

#endif
