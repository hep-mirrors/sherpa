#ifndef AMISIC_Tools_NonPerturbative_XSecs_H
#define AMISIC_Tools_NonPerturbative_XSecs_H

#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"


namespace REMNANTS {
  class Remnant_Handler;
}

namespace AMISIC {
  struct event_mode {
    enum code {
      elastic = 1,
      SDA     = 2,
      SDB     = 3,
      DD      = 4,
      unknown = 99
    };
  };

  
  class NonPerturbative_XSecs {
  private:
    Hadronic_XSec_Calculator    * p_xsecs;
    evt_type::code                m_evttype;
    std::array<ATOOLS::Flavour,2> m_inflav;
    std::array<ATOOLS::Flavour,4> m_outflav;
    double m_masses[2], m_masses2[2], m_outmasses[4], m_outmasses2[4];
    double m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_alphaQED;
    double m_s, m_s0, m_mpi, m_deltaMres, m_cres;
    ATOOLS::Vec4D m_inmom[2];
    ATOOLS::Vec4D m_outmom[4];

    inline const double beta(const size_t hadtag,const double & t) const {
      return sqrt(p_xsecs->s_X[hadtag][hadtag]) * exp(p_xsecs->s_slopes[hadtag]*t);
    }
    inline const double F_SD(const double & M2,const double & Mres2=4.,const double & cres=2.) const {
      return ( (1.-M2/m_s) * (1.+cres*Mres2/(Mres2+M2)) );
    }
    inline const double ExponentialDist(const double & smin,const double & smax,
					const double & pref) {
      double ran = ATOOLS::ran->Get(); 
      return log(ran*exp(pref*smax) + (1.-ran)*exp(pref*smin))/pref;
    }
    inline const double ExponentialWeight(const double & smin,const double & smax,const double & pref) {
      return (exp(pref*smax) - exp(pref*smin)) / pref;
    }
    inline const double PowerDist(const double & smin,const double & smax,
				  const double & expo) {
      double ran = ATOOLS::ran->Get(); 
      if (expo==-1.) return smin*pow(smax/smin,ran);
      return pow(ran*pow(smax,expo+1.)+(1.-ran)*pow(smin,expo+1.), 1./(expo+1.));
    }
    inline const double PowerWeight(const double & smin,const double & smax,const double & expo) {
      if (expo==-1.) return log(smax/smin);
      return (pow(smax,expo+1.)-pow(smin,expo+1.))/(expo+1.);
    }
    
    const event_mode::code SelectMode(std::array<size_t,2> & hadtags);
    double DiffElXSec(const std::array<size_t,2> hadtags,const double & s,const double & t) const;
    bool   SplitDiffractiveState(const size_t & pos);
    ATOOLS::Blob * ElasticScatter(REMNANTS::Remnant_Handler * remnants,
				  const std::array<size_t,2> hadtags);
    ATOOLS::Blob * SingleDiffractiveScatter(REMNANTS::Remnant_Handler * remnants,
					    const std::array<size_t,2> hadtags,const size_t & diffbeam);
    ATOOLS::Blob * DoubleDiffractiveScatter(REMNANTS::Remnant_Handler * remnants,
					    const std::array<size_t,2> hadtags);
    ATOOLS::Blob * InitBlob(REMNANTS::Remnant_Handler * remnants,
			    const double & muR=0.,const double & muQ=0.);
  public:
    NonPerturbative_XSecs(Hadronic_XSec_Calculator * xsecs);

    ATOOLS::Blob * MakeScatter(REMNANTS::Remnant_Handler * remnants);
  };
}

#endif
