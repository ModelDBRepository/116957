/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_BASE_H
#define CN_BASE_H

#include <cassert>
#include <iostream>
#include <list>
#include <vector>

#define forall(l, it) for (it= (l).begin(); it != (l).end(); it++)
#define pw2(x) x*x
#define pw3(x) x*x*x
#define pw4(x) x*x*x*x

#define NEURTYPENO 33

#define TIMENEURON -1
#define HHNEURON 0
#define HHCANEURON 1
#define FNGMNEURON 2
#define ICANEURON 3
#define STUPIDNEURON 4
#define POISSONNEURON 5
#define VALNEURON 6
#define MULTIFIRE_INPUTNEURON 7
#define IFNEURON 8
#define KOLNEURON 9
#define KOLINNEURON 10
#define KOLMULTIFIRE_INPUTNEURON 11
#define LPNEURON 12
#define LPGNEURON 13
#define LPJNEURON 14
#define LPTNEURON 15
#define HVCE1 16
#define HVCI1 17
#define ECNEURON 18
#define LPANEURON 19
#define LPMNEURON 20
#define HRNEURON 21
#define LPRNEURON 22
#define PNRAMON 23
#define LNRAMON 24
#define PSEUDONEURON 25
#define KCDNEURON 26
#define LMPNEURON 27
#define LTVNEURON 28
#define VDPOLNEURON 29
#define SIN 30
#define DATA 31
#define PNANEURON 32

#define SYNTYPENO 31

#define DCINPUT 0
#define DEMIGAP 1
#define RALL 2
#define ALINSYN 3
#define CRALL 4
#define DYNSTDP 5
#define HERA 6
#define LRNRALL 7
#define GRAD 8
#define KOLSYNAPSE 9
#define HEBBKOL 10
#define KOLGRADSYNAPSE 11
#define HVCSYN 12
#define ABSYN 13
#define ABECPLAST 14
#define SYNAS 15
#define ABECPLAST3 16
#define SYNASPLAST 17
#define RMSYN 18
#define sRMSYN 19
#define IRMSYN 20
#define SMSTDP1 21
#define LTVSYN 22
#define INPUTFUNCTION 23
#define S01SYN 24
#define S01ECPLAST3 25
#define RALLECPLAST3 26
#define T2RALL 27
#define T2RALLECPLAST3 28
#define ABECHEBB3 29
#define DEMIGAPSYNAPSE 30

#endif
