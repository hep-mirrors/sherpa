#ifndef SHRIMPS_Eikonals_Rapidity_Density_H
#define SHRIMPS_Eikonals_Rapidity_Density_H

#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"

namespace SHRIMPS {
  class Rapidity_Density: public ATOOLS::Function_Base {
  private:
    double m_Delta, m_lambda;
    double m_b1, m_b2, m_max, m_mean;
    Eikonal_Contributor * p_omegaik, * p_omegaki;

    double Integrate(const double & ymin,const double & ymax);
    double MaxWeight();
    double DeltaOmega(const double & y1,const double & y2);
  public:
    Rapidity_Density(const double & Delta,const double & lambda);
    ~Rapidity_Density() {}
    void   SetEikonal(Omega_ik * eikonal);
    void   SetImpactParameters(const double & b1, const double & b2);
    double operator()(double y);
    size_t NGluons(const double & ymin,const double & ymax);
    double SelectRapidity(const double & ymin,const double & ymax);
    double SingletWeight(const double & y1,const double & y2);
    double OctetWeight(const double & y1,const double & y2);

    const double MeanNGluons() const { return m_mean; }
    void  SetLambdaForTest(const double & lambda) { m_lambda = lambda; } 
  };
}
#endif
