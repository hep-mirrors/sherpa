#ifndef REMNANTS_Tools_Primordial_KPerp_H
#define REMNANTS_Tools_Primordial_KPerp_H

#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Phys/Particle.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"
#include <array>
#include <map>
#include <string>

namespace REMNANTS {
  class Remnant_Handler;

  class Primordial_KPerp {
  private:
    ATOOLS::Flavour m_beamflav;
    std::array<primkT_form, 2>   m_form;
    std::array<primkT_recoil, 2> m_recoil;
    std::array<double, 2> m_SIMean, m_SISigma, m_SIQ2, m_SpecMean, m_SpecSigma,
            m_SpecQ2;
    std::array<double, 2>   m_SIKtmax, m_SIEta, m_SpecKtmax, m_SpecEta;
    double                  m_mean, m_sigma, m_Q2, m_ktmax, m_eta;
    size_t                  m_beam;
    ATOOLS::ParticleMomMap* p_ktmap;

    double KT_Gauss(const double & ktmax) const;
    double KT_Gauss_Limited(const double & ktmax) const;
    double KT_Dipole(const double & ktmax) const;
    double KT_Dipole_Limited(const double & ktmax) const;
    void   BalanceKT(const ATOOLS::Vec4D & kt_Show,const double & E_Show,
		     const ATOOLS::Vec4D & kt_Spec,const double & E_Spec);
    double DipoleWeight(const double & kt) const;
    double LimitedWeight(const double & kt) const;

    bool                                      m_analysis;
    std::map<std::string, ATOOLS::Histogram * >m_histos;
    void InitAnalysis();
    void FinishAnalysis();
  public:
    Primordial_KPerp();
    ~Primordial_KPerp();

    void Initialize(Remnant_Handler * rhandler);
    bool CreateBreakupKinematics(const size_t & beam,ATOOLS::ParticleMomMap * ktmap,
				 const double & scale);
    ATOOLS::Vec4D KT(const ATOOLS::Particle * part,const double & ktext=-1.);
  };
}

#endif
