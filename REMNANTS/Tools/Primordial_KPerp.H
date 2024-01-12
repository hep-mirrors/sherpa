#ifndef REMNANTS_Tools_Primordial_KPerp_H
#define REMNANTS_Tools_Primordial_KPerp_H

#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace REMNANTS {
  class Remnant_Handler;
  
  class Primordial_KPerp {
  private:
    pkform::code   m_form[2];
    pkrecoil::code m_recoil[2];
    double m_SIMean[2], m_SISigma[2], m_SIQ2[2], m_SpecMean[2], m_SpecSigma[2], m_SpecQ2[2];
    double m_SIKtmax[2], m_SIEta[2], m_SpecKtmax[2], m_SpecEta[2];
    double m_mean, m_sigma, m_Q2, m_ktmax, m_eta;
    size_t m_beam;
    ATOOLS::ParticleMomMap * p_ktmap;

    double KT_Gauss(const double & ktmax) const;
    double KT_Gauss_Limited(const double & ktmax) const;
    double KT_Dipole(const double & ktmax) const;
    double KT_Dipole_Limited(const double & ktmax) const;
    void   BalanceKT(const ATOOLS::Vec4D & kt_Show,const double & E_Show,
		     const ATOOLS::Vec4D & kt_Spec,const double & E_Spec);
    double DipoleWeight(const double & kt) const;
    double LimitedWeight(const double & kt) const;

    bool m_analysis;    
    std::map<std::string, ATOOLS::Histogram * >m_histos;
    void InitAnalysis();
    void FinishAnalysis();
  public:
    Primordial_KPerp();
    ~Primordial_KPerp();

    void Initialize(Remnant_Handler * rhandler);
    bool CreateBreakupKinematics(const size_t & beam,ATOOLS::ParticleMomMap * ktmap,
				 const double & scale);
    ATOOLS::Vec4D KT(const ATOOLS::Particle * part);
  };
}

#endif
