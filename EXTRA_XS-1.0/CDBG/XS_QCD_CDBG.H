#ifndef XS_QCD_CDBG_H
#define XS_QCD_CDBG_H

#include "Single_XS.H"

namespace EXTRAXS {

  class N_Gluon_CDBG;

  class CDXS_gg_ng: public Single_XS {
  private:

    N_Gluon_CDBG *p_bg;

    std::vector<ATOOLS::Vec4D> m_p;
  
  public:

    // constructor
    CDXS_gg_ng(const size_t &nin,const size_t &nout,
	       const ATOOLS::Flavour *const fl);
    
    // destructor
    ~CDXS_gg_ng();

    double operator()(const double s,const double t,const double u);
    bool   SetColours(const double s,const double t,const double u);

    bool Tests();

  };
 
}// end of namespace EXTRAXS

#endif

