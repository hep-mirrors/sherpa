#ifndef SINIC_Cross_Sections_Cross_Sections_H
#define SINIC_Cross_Sections_Cross_Sections_H

#include "SINIC++/Cross_Sections/Sigma_Inelastic.H"
#include "SINIC++/Cross_Sections/Sigma_Elastic.H"
#include "SINIC++/Cross_Sections/Sigma_SD.H"
#include "SINIC++/Cross_Sections/Sigma_DD.H"
#include "SINIC++/Cross_Sections/Sigma_Total.H"
#include "SINIC++/Eikonals/Omega_ik.H"
#include <list>
#include <map>

namespace SINIC {
  class Cross_Sections {
  private:
    std::list<Omega_ik *> * p_eikonals;    
    Omega_ik              * p_selected;       
    double                  m_originalY, m_cutoffY, m_Y; 
    double                  m_B, m_Bmin, m_Bmax;
    Sigma_Inelastic         m_sigma_inelastic;
    Sigma_Elastic           m_sigma_elastic;
    Sigma_SD                m_sigma_SD;
    Sigma_DD                m_sigma_DD;
    double                  m_stot, m_sinel, m_sel, m_sSD, m_sDD, m_sval;
    std::map<int,double>    m_modemap;
    int                     m_test;
    double 	            m_Delta,m_pref,m_Lambda2,m_beta,m_kappa;
  public:
    Cross_Sections() {}
    Cross_Sections(std::list<Omega_ik *> * eikonals,const double & energy,
		   const int & test);
    ~Cross_Sections();
    
    void CalculateTotalCrossSections();
    int  SelectCollisionMode();

    std::list<Omega_ik *> * GetEikonals() { return p_eikonals; }    

    Sigma_Elastic   * GetSigmaElastic()   { return (&m_sigma_elastic); }
    Sigma_SD        * GetSigmaSD()        { return (&m_sigma_SD); }
    Sigma_DD        * GetSigmaDD()        { return (&m_sigma_DD); }
    Sigma_Inelastic * GetSigmaInelastic() { return (&m_sigma_inelastic); }

    const double & SigmaTot() const     { return m_stot; }
    const double & SigmaSD() const      { return m_sSD; }
    const double & SigmaDD() const      { return m_sDD; }
    const double & SigmaEl() const      { return m_sel; }
    const double & SigmaInel() const    { return m_sinel; }
    const double & ElasticSlope() const { return m_sval; }
  };
}
#endif
