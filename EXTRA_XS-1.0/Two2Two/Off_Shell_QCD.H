#ifndef Off_Shell_QCD_H
#define Off_Shell_QCD_H

#include "Single_XS.H"

namespace EXTRAXS {

  class Off_Shell_gg_qqb: public Single_XS {  
  private:
    
    double m_alphas;

  public:

    // constructor
    Off_Shell_gg_qqb(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		     const int scalescheme,const int kfactorscheme,const double scalefactor);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

  };// end of class Off_Shell_gg_qqb
 
}// end of namespace EXTRAXS

#endif
