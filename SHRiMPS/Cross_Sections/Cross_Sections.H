#ifndef SHRIMPS_Cross_Sections_Cross_Sections_H
#define SHRIMPS_Cross_Sections_Cross_Sections_H

#include "SHRiMPS/Eikonals/Omega_ik.H"
#include <list>
#include <map>

namespace SHRIMPS {
  class Cross_Sections {
  private:
    Omega_ik              * p_selected;       
    double m_xstot, m_slope, m_xsinel, m_xsel, m_xsSD, m_xsDD;
    //Sigma_Inelastic m_sigma_inelastic;
    //Sigma_Elastic   m_sigma_elastic;
    //Sigma_SD        m_sigma_SD;
    //Sigma_DD        m_sigma_DD;
    std::map<run_mode::code,double> m_modemap;
  public:
    Cross_Sections();
    ~Cross_Sections();
    
    void CalculateCrossSections();
    void Test(const std::string & dirname);
    //run_mode::code SelectCollisionMode();

    //std::list<Omega_ik *> * GetEikonals() { return p_eikonals; }    

    //Sigma_Elastic   * GetSigmaElastic()   { return (&m_sigma_elastic); }
    //Sigma_SD        * GetSigmaSD()        { return (&m_sigma_SD); }
    //Sigma_DD        * GetSigmaDD()        { return (&m_sigma_DD); }
    //Sigma_Inelastic * GetSigmaInelastic() { return (&m_sigma_inelastic); }

    const double & SigmaTot() const     { return m_xstot; }
    const double & SigmaSD() const      { return m_xsSD; }
    const double & SigmaDD() const      { return m_xsDD; }
    const double & SigmaEl() const      { return m_xsel; }
    const double & SigmaInel() const    { return m_xsinel; }
    const double & ElasticSlope() const { return m_slope; }
  };
}
#endif
