#ifndef PHASIC_Process_Spin_Color_Correlated_ME2_H
#define PHASIC_Process_Spin_Color_Correlated_ME2_H

#include "PHASIC++/Process/External_ME_Args.H"
#include "MODEL/Main/Model_Base.H"

namespace MODEL { struct Coupling_Map; }

namespace PHASIC {

  /* Spin- and Color-correlated matrix elements [arXiv:hep-ph/9605323]
     of the form <1,2,3,... | Vij*Ti*Tj |...,3,2,1 > */

  class Spin_Color_Correlated_ME2 {

  protected:

    MODEL::Coupling_Data* p_aqcd, * p_aqed;

  public:

    Spin_Color_Correlated_ME2(const External_ME_Args& args);

    virtual ~Spin_Color_Correlated_ME2() {}

    /* Calculate <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
    virtual double CalcColorCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
				       const size_t& born_ij,
				       const size_t& born_k) const = 0;

    /* Calculate <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
    virtual double CalcSpinCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
				      const ATOOLS::Vec4D& p_tilde,
				      const size_t& born_ij,
				      const size_t& born_k) const = 0;

    void SetCouplings(MODEL::Coupling_Data* p_rqcd,
		      MODEL::Coupling_Data* p_rqed);
    
    virtual double AlphaQCD() const;
    virtual double AlphaQED() const;

    /* Need access methods in order to synchronize these couplings
       with wrapper processes */
    MODEL::Coupling_Data* RunningQCD() const { return p_aqcd; }
    MODEL::Coupling_Data* RunningQED() const { return p_aqed; }

    static Spin_Color_Correlated_ME2 *GetME2(const External_ME_Args& args);

  };

}


#define DECLARE_SPINCOLORCORRELATEDME2_GETTER(NAME,TAG)		                         \
  DECLARE_GETTER(NAME,TAG,PHASIC::Spin_Color_Correlated_ME2,PHASIC::External_ME_Args);   \
  void ATOOLS::Getter<PHASIC::Spin_Color_Correlated_ME2,PHASIC::External_ME_Args,NAME>:: \
  PrintInfo(std::ostream &str,const size_t width) const		                         \
  {                                                                                      \
    str<<#TAG;                                                                           \
  }

#endif
