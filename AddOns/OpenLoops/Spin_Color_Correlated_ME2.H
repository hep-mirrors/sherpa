#ifndef AddOns_OpenLoops_Spin_Color_Correlated_ME2_H
#define AddOns_OpenLoops_Spin_Color_Correlated_ME2_H

#include "PHASIC++/Process/Spin_Color_Correlated_ME2.H"

namespace OpenLoops {
  
  class Spin_Color_Correlated_ME2 : public PHASIC::Spin_Color_Correlated_ME2 {

  public:
    
    Spin_Color_Correlated_ME2(const PHASIC::Correlator_Args& args,
			      int ol_id, const AmplitudeType& type) :
      PHASIC::Spin_Color_Correlated_ME2(args),
      m_ol_id(ol_id), m_amptype(type) {}


    ~Spin_Color_Correlated_ME2() {}

    double CalcColorCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
			       const size_t& born_ij,
			       const size_t& born_k) const;

    double CalcSpinCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
			      const ATOOLS::Vec4D& p_tilde,
			      const size_t& born_ij,
			      const size_t& born_k) const;

  private:

    int m_ol_id;
    AmplitudeType m_amptype;
    const AmplitudeType AmpType() const { return m_amptype; }

  };

}

#endif
