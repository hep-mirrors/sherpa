#ifndef SHRIMPS_Cross_Sections_Sigma_Base_H
#define SHRIMPS_Cross_Sections_Sigma_Base_H
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Sigma_Base : public ATOOLS::Function_Base {
  protected:
    std::list<Omega_ik *> * p_eikonals;
    Omega_ik              * p_eikonal;
    double                  m_originalY, m_cutoffY, m_Y, m_E;
    double                  m_accu;
    double                  m_sigma;
  public:
    Sigma_Base() : p_eikonals(MBpars.GetEikonals()), m_sigma(0.) {}

    virtual void   SetEikonal(Omega_ik * eikonal)    { p_eikonal = eikonal; }
    virtual void   SetSigma(const double & sigma=-1) { m_sigma = sigma; }
    virtual const double & XSec() const              { return m_sigma; }

    virtual double Calculate(Omega_ik * eikonal=NULL);
    virtual double operator()(double B);

    virtual double GetValue(const double & B)         = 0;
    virtual double GetCombinedValue(const double & B) = 0;
  };
}
#endif
