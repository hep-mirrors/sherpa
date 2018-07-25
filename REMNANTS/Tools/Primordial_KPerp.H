#ifndef REMNANTS_Tools_Primordial_KPerp_H
#define REMNANTS_Tools_Primordial_KPerp_H

#include "REMNANTS/Tools/Remnant_Parameters.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace REMNANTS {
  class Primordial_KPerp {
  private:
    prim_kperp_form::code m_form[2];
    double m_mean[2], m_sigma[2], m_Q2[2], m_ktmax[2], m_eta[2];
    size_t m_beam;
    ATOOLS::ParticleMomMap * p_ktmap;

    prim_kperp_form::code SelectForm(const std::string & form=std::string("dipole_limited"));
    ATOOLS::Vec4D         KT(const ATOOLS::Vec4D & mom);
    double KT_Gauss(const double & ktmax) const;
    double KT_Gauss_Limited(const double & ktmax) const;
    double KT_Dipole(const double & ktmax) const;
    double KT_Dipole_Limited(const double & ktmax) const;
    void   BalanceKT(const ATOOLS::Vec4D & kt_tot,const double & E_tot);
    double DipoleWeight(const double & kt) const;
    double LimitedWeight(const double & kt) const;

    bool m_analysis;    
    std::map<std::string, ATOOLS::Histogram * >m_histos;
    void InitAnalysis();
    void FinishAnalysis();
  public:
    Primordial_KPerp();
    ~Primordial_KPerp();

    void Initialize();
    bool CreateBreakupKinematics(const size_t & beam,ATOOLS::ParticleMomMap * ktmap,
				 const double & scale);
  };
}

#endif
