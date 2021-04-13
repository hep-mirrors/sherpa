#ifndef SHRIMPS_Eikonals_Rapidity_Density_H
#define SHRIMPS_Eikonals_Rapidity_Density_H

#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"

namespace SHRIMPS {
  class Rapidity_Density: public ATOOLS::Function_Base {
  private:
    double              m_Delta, m_lambda, m_Ymax;
    absorption::code    m_absorp;
    double              m_b1, m_b2, m_max, m_mean;
    Eikonal_Contributor * p_omegaik, * p_omegaki;

    double Integrate(const double & ymin,const double & ymax);
    double DeltaOmega(const double & y1,const double & y2);
  public:
    Rapidity_Density(const double & Delta,const double & lambda,const double & Ymax,
		     const absorption::code & absorp);
    ~Rapidity_Density() {}
    void   SetEikonal(Omega_ik * eikonal);
    void   SetImpactParameters(const double & b1, const double & b2);
    double operator()(double y);
    size_t NGluons(const double & ymin,const double & ymax);
    double SelectRapidity(const double & ymin,const double & ymax);
    double AbsorptionWeight(double y);
    double SingletWeight(const double & y1,const double & y2);
    double OctetWeight(const double & y1,const double & y2);
    double EffectiveIntercept(const double & b1, const double & b2,const double & y);
    double MaxWeight();
    
    const double MeanNGluons() const              { return m_mean; }
    void  SetLambdaForTest(const double & lambda) { m_lambda = lambda; } 
    void  Test(Omega_ik * eikonal);
  };
}
#endif
