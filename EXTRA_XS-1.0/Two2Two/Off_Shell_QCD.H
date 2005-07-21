#ifndef Off_Shell_QCD_H
#define Off_Shell_QCD_H

#include "Single_XS.H"
#include "Info_Key.H"

namespace EXTRAXS {

  class Off_Shell_gg_qqb: public Single_XS {  
  private:
    
    double m_alphas;

  public:

    // constructor
    Off_Shell_gg_qqb(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

  };// end of class Off_Shell_gg_qqb
 
  class Off_Shell_gg_gg: public Single_XS {  
  private:
    
    double m_alphas, m_qtcut;

    ATOOLS::Info_Key m_zkey[2];

    unsigned int m_jets;

  public:

    // constructor
    Off_Shell_gg_gg(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    bool Trigger(const ATOOLS::Vec4D *const momenta);

  };// end of class Off_Shell_gg_gg
 
}// end of namespace EXTRAXS

#endif
