#ifndef XS_QCD_CDBG_T_H
#define XS_QCD_CDBG_T_H

#include "Single_XS.H"

namespace EXTRAXS {

  class N_Parton_CDBG_T;

  class CDXS_T_pp_np: public Single_XS {
  private:

    N_Parton_CDBG_T *p_bg;

    std::vector<ATOOLS::Vec4D> m_p;
  
  public:

    // constructor
    CDXS_T_pp_np(const size_t &nin,const size_t &nout,
	       const ATOOLS::Flavour *const fl);
    
    // destructor
    ~CDXS_T_pp_np();

    double operator()(const double s,const double t,const double u);
    bool   SetColours(const double s,const double t,const double u);

    bool Tests();

  };
 
}// end of namespace EXTRAXS

#endif

