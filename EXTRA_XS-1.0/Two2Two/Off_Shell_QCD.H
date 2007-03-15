#ifndef Off_Shell_QCD_H
#define Off_Shell_QCD_H

#include "Single_XS.H"
#include "Info_Key.H"

namespace PDF { class Doubly_Unintegrated_PDF; }

namespace EXTRAXS {

  class Off_Shell_gg_qqb: public Single_XS {  
  private:
    
    double m_alphas;

    PDF::Doubly_Unintegrated_PDF *p_pdfs[2];

  public:

    // constructor
    Off_Shell_gg_qqb(const size_t nin,const size_t nout,
		     const ATOOLS::Flavour *flavours);
    
    // member functions
    double operator()(double s,double t,double u);
    bool   SetColours(double s,double t,double u);

    bool SetPDFFlavours();

  };// end of class Off_Shell_gg_qqb
 
}// end of namespace EXTRAXS

#endif
