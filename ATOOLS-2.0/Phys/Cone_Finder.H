#ifndef Cone_Finder_h
#define Cone_Finder_h

#include "Parton_List.H"
#include "Selector.H"
#include "Poincare.H"

namespace APHYTOOLS {
  class Cone_Finder : public Selector_Base {
    double m_rcone;
    double * m_value;
        
    double Rmin(AMATOOLS::Vec4D *);

    double DEta12(AMATOOLS::Vec4D &,AMATOOLS::Vec4D &);
    double DPhi12(AMATOOLS::Vec4D &,AMATOOLS::Vec4D &);
  public:
    Cone_Finder(int,double);
    
    void   Init(const AMATOOLS::Vec4D *);
    bool   Trigger(const AMATOOLS::Vec4D *);
    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);
    double * ActualValue();

  };
}

#endif








