#ifndef Amplitude_Manipulator_h
#define Amplitude_Manipulator_h

#include "Flavour.H"
#include "Topology.H"
#include "Single_Amplitude.H"

namespace AMEGIC {

  class Amplitude_Manipulator {
    int N;
    APHYTOOLS::Flavour* fl;
    int* b;
    
    void   SetPrev(Point*);
    int    Permutation(int*,int);
    void   GetPerm(int*,Single_Amplitude*,int&);
    Point* FindNext(Point*);
    void   GetFermionLine(Point*,Point*&,Point*&);
    Point* ForwardLine(Point*);
    Point* BackwardLine(Point*);
    
    int  SetPropOrientation(Point*,Point*);
    void ForwardLineOrientation(Point*,int&);
    void BackwardLineOrientation(Point*,int&);
    
    void SetFermionNumberFlow(Point*,Point*);
    void SetForwardFNFlow(Point*,int);
    void SetBackwardFNFlow(Point*,int);
    
  public:
    Amplitude_Manipulator(int _no,APHYTOOLS::Flavour* _fl,int* _b): N(_no), fl(_fl), b(_b) {}
    void   FixSign(Single_Amplitude*);  
  };
}
#endif











