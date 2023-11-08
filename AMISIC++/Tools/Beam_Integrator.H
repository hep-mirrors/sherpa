#ifndef AMISIC_Tools_Beam_Integrator_H
#define AMISIC_Tools_Beam_Integrator_H

#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Random.H"

namespace AMISIC {
  class Beam_Integrator {
  private:
    Hadronic_XSec_Calculator   * p_xsecs;
    BEAM::Beam_Base            * p_beams[2];
    std::array<ATOOLS::Vec4D, 2> m_beammoms,m_inmoms;
    
    size_t         m_Npoints, m_variable;
    int            m_sign;
    double         m_smin, m_smax, m_swt;
    double         m_sprime, m_yprime, m_total, m_error, m_max;
    evt_type::code m_evttype;

    void   CalculateXSecs();
    double MakePoint();
    double TrialEvent();
    double inline MEweight() {
      double MEweight     = 0.;
      if (m_evttype==evt_type::Elastic       ||
	  m_evttype==evt_type::QuasiElastic) MEweight += p_xsecs->SigmaEl();
      if (m_evttype==evt_type::DiffractiveA  ||
	  m_evttype==evt_type::QuasiElastic) MEweight += p_xsecs->SigmaSDA();
      if (m_evttype==evt_type::DiffractiveB  ||
	  m_evttype==evt_type::QuasiElastic) MEweight += p_xsecs->SigmaSDB();
      if (m_evttype==evt_type::DiffractiveAB ||
	  m_evttype==evt_type::QuasiElastic) MEweight += p_xsecs->SigmaDD();
      return MEweight;
    }
  public:
    Beam_Integrator() : m_Npoints(10000), m_variable(0), m_smin(4.),
			m_total(0.), m_error(1.), m_max(0.) {}
    void Init(Hadronic_XSec_Calculator * xsecs,
	      BEAM::Beam_Base * beam1,BEAM::Beam_Base * beam2,
	      evt_type::code & type);
    bool operator()();
    inline const double        & Sprime()                    const { return m_sprime; }
    inline const double        & Yprime()                    const { return m_yprime; }
    inline const ATOOLS::Vec4D & Momentum(const size_t beam) const { return m_inmoms[beam]; }
  };
}

#endif
