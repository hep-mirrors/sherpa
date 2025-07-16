#ifndef MODEL_Main_Effective_Higgs_Coupling_H
#define MODEL_Main_Effective_Higgs_Coupling_H

#include "ATOOLS/Math/MyComplex.H"

namespace MODEL {


  class Effective_Higgs_Coupling  {
    double m_mhiggs;
    Complex f(double); 
  public:

    // constructor
    Effective_Higgs_Coupling(const double mhiggs);
    
    Complex GetFermionContribution(double mass,bool pseudoscalar=false);
    Complex GetVectorContribution(double mass);
    Complex GetScalarContribution(double mass);
  };

}



#endif
