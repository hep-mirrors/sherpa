#ifndef Cone_Finder_h
#define Cone_Finder_h

#include "Parton_List.H"
#include "Selector.H"
#include "Poincare.H"

namespace ATOOLS {
  class Cone_Finder : public Selector_Base {
    double m_rcone;
    double * m_value;
        
    double Rmin(Vec4D *);

    double DEta12(Vec4D &,Vec4D &);
    double DPhi12(Vec4D &,Vec4D &);
  public:
    Cone_Finder(int,Flavour*,double);
    
    void   Init(const Vec4D *);
    bool   Trigger(const Vec4D *);
    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();

  };
}

#endif








