#ifndef running_fermion_mass_h
#define running_fermion_mass_h

#include "Flavour.H"
#include "Running_AlphaS.H"
#include "Function_Base.H"

namespace MODEL {
  class Running_Fermion_Mass : public AMATOOLS::Function_Base {
  protected:
    int              m_order;
    double           m_gamma,m_beta;
    double           m_polemass,m_alphahat;
    Running_AlphaS * p_as;
  public:
    Running_Fermion_Mass(APHYTOOLS::Flavour,double,Running_AlphaS *);
    ~Running_Fermion_Mass() {}
    double operator()(double); 
    void SelfTest();
  };
}



#endif
