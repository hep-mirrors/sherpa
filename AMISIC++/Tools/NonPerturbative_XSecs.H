#ifndef AMISIC_Tools_NonPerturbative_XSecs_H
#define AMISIC_Tools_NonPerturbative_XSecs_H

#include "AMISIC++/Tools/Beam_Integrator.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Random.H"


namespace REMNANTS {
  class Remnant_Handler;
}

namespace BEAM {
  class Beam_Base;
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
    REMNANTS::Remnant_Handler   * p_remnants;
    Hadronic_XSec_Calculator    * p_xsecs;
    BEAM::Beam_Base             * p_beams[2];
    bool                          m_variableS;
    Beam_Integrator               m_integrator;
    ATOOLS::Poincare              m_boost;

    evt_type::code                m_evttype;
    std::array<ATOOLS::Flavour,2> m_inflav;
    std::array<double,2>          m_masses, m_masses2;
    std::array<ATOOLS::Flavour,4> m_outflav;
    std::array<double,4>          m_outmasses, m_outmasses2;
    std::array<size_t,2>          m_hadtags;
    ATOOLS::Vec4D                 m_inmom[2], m_outmom[4];
    double m_smin, m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_alphaQED;
    double m_s, m_s0, m_mpi, m_mpi2, m_deltaMres, m_cres;
    double m_mrho, m_mrho2, m_mrho_min, m_momega, m_momega2;
    double m_Grho, m_Grho2, m_Gomega, m_Gomega2;
    double m_q_rho, m_q_omega, m_f_omega, m_phi_omega;
    double m_f_nr, m_Lambda2_nr, m_delta_nr, m_A2max;
    double m_weight, m_xsec, m_muf2, m_mur2;
    two_pions::code m_twopions;

    unsigned long int m_calls, m_fails;

    inline double beta(const size_t hadtag,const double & t) const {
      return sqrt(Hadronic_XSec_Calculator::s_X[hadtag][hadtag]) *
             exp(Hadronic_XSec_Calculator::s_slopes[hadtag]*t);
    }
    inline double F_SD(const double & M2,const double & Mres2) const {
      return ( (1.-M2/m_s) * (1.+m_cres*Mres2/(Mres2+M2)) );
    }
    inline double F_DD(const std::array<double, 2> & M2,
			     const std::array<double, 2> & Mres2) const {
      double averMass = sqrt(M2[0])+sqrt(M2[1])/2., Mass2 = 4.*averMass*averMass;
      return ( (1.-Mass2/m_s) *
	       (m_masses[0]*m_masses[1]*m_s)/(m_masses[0]*m_masses[1]*m_s+M2[0]*M2[1]) *
	       (1.+m_cres*Mres2[0]/(Mres2[0]+M2[0])) *
	       (1.+m_cres*Mres2[1]/(Mres2[1]+M2[1])) );
    }
    inline ATOOLS::Flavour GetFlavour(const size_t & pos,const size_t & had) const {
      if (m_inflav[pos]==ATOOLS::Flavour(kf_p_plus) ||
	  m_inflav[pos]==ATOOLS::Flavour(kf_p_plus).Bar() ||
	  m_inflav[pos]==ATOOLS::Flavour(kf_n) ||
	  m_inflav[pos]==ATOOLS::Flavour(kf_n).Bar()) return m_inflav[pos];
      if (m_inflav[pos]==ATOOLS::Flavour(kf_photon)) {
	switch (had) {
	case 0:  return ATOOLS::Flavour(kf_rho_770);
	case 1:  return ATOOLS::Flavour(kf_omega_782);
	case 2:  return ATOOLS::Flavour(kf_phi_1020);
	case 3:  return ATOOLS::Flavour(kf_J_psi_1S);
	default: break;
	}
      }
      return ATOOLS::Flavour(kf_none);
    }

    inline const double ExponentialDist(const double & smin,const double & smax,
					const double & pref) {
      double ran = ATOOLS::ran->Get();
      return log(ran*exp(pref*smax) + (1.-ran)*exp(pref*smin))/pref;
    }
    inline const double ExponentialWeight(const double & smin,const double & smax,
					  const double & pref) {
      return (exp(pref*smax) - exp(pref*smin)) / pref;
    }
    inline const double PowerDist(const double & smin,const double & smax,
				  const double & expo) {
      double ran = ATOOLS::ran->Get();
      if (expo==-1.) return smin*pow(smax/smin,ran);
      return pow(ran*pow(smax,expo+1.)+(1.-ran)*pow(smin,expo+1.), 1./(expo+1.));
    }
    inline const double PowerWeight(const double & smin,const double & smax,
				    const double & expo) {
      if (expo==-1.) return log(smax/smin);
      return (pow(smax,expo+1.)-pow(smin,expo+1.))/(expo+1.);
    }

    const event_mode::code SelectMode();
    bool   FixFS(std::array<ATOOLS::Flavour, 2> & flavs);
    bool   FixOutMomenta(const double & t);
    bool   SplitDiffractiveState(const size_t & pos);
    bool   SplitDiffractiveStateIntoHadrons(const size_t & pos,const double & M2);
    bool   SelectFlavoursOfDiffraction(const size_t & pos,const double & M2);
    bool   SetRhoMasses2();
    double RhoMassModifier(const double & M2);
    void   SplitRhoIntoPions(const size_t & pos);
    bool   ElasticScatter(ATOOLS::Blob * blob);
    bool   SingleDiffractiveScatter(const size_t & pos,ATOOLS::Blob * blob);
    bool   DoubleDiffractiveScatter(ATOOLS::Blob * blob);
    void   InitBlob(ATOOLS::Blob * blob,const double & muR2=0.,const double & muQ2=0.);

    bool   m_ana;
    size_t m_Nsteps, m_Nevents;
    double m_maxT, m_minM, m_maxM;
    std::map<std::string,ATOOLS::Histogram *> m_histos;

    void InitHistos();
    void WriteHistos();
    void Tests();
    void TestElastic();
    void TestSingleD();
    void TestDoubleD();
  public:
    NonPerturbative_XSecs(REMNANTS::Remnant_Handler * remnants,
			  Hadronic_XSec_Calculator * xsecs);
    ~NonPerturbative_XSecs();
    void SetBeams(BEAM::Beam_Base * beam1,BEAM::Beam_Base * beam2);
    void CalculateSDependentCrossSections();

    bool MakeScatter(ATOOLS::Blob * blob);

    inline const double & Weight() const { return m_weight; }
    inline const double & Smin()   const { return m_smin; }
    inline const double & MuF2()   const { return m_muf2; }
    inline const double & MuR2()   const { return m_mur2; }

    inline void SetHadTags(const std::array<size_t,2> tags) {
      for (size_t i=0;i<2;i++) m_hadtags[i] = tags[i];
    }
    double DiffElXSec(const double & s,const double & t) const;
    double DiffSDXSec(const double & s,const double & t,const double & M2,
		      const size_t & pos) const;
    double DiffDDXSec(const double & s,const double & t,
		      const std::array<double, 2> & M2) const;
  };
}

#endif
