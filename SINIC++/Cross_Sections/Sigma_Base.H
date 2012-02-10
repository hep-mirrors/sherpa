#ifndef SINIC_Cross_Sections_Sigma_Base_H
#define SINIC_Cross_Sections_Sigma_Base_H
#include "SINIC++/Eikonals/Omega_ik.H"
#include "SINIC++/Tools/MinBias_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SINIC {
  class Sigma_Base : public ATOOLS::Function_Base {
  protected:
    std::list<Omega_ik *> * p_eikonals;
    Omega_ik              * p_eikonal;
    double                  m_originalY, m_cutoffY, m_Y, m_E;
    double                  m_accu;
    double                  m_sigma;
  public:
    Sigma_Base() {}
    Sigma_Base(std::list<Omega_ik *> * eikonals) : 
      p_eikonals(eikonals), 
      m_originalY(MBpars("originalY")), 
      m_cutoffY(MBpars("deltaY")), m_Y(m_originalY-m_cutoffY),
      m_E(ATOOLS::Flavour(kf_p_plus).HadMass()*exp(m_originalY)),
      m_accu(MBpars("accu")), m_sigma(0.) {}
    virtual void   SetEikonal(Omega_ik * eikonal)    { p_eikonal = eikonal; }
    virtual void   SetSigma(const double & sigma=-1) { m_sigma = sigma; }
    virtual double Calculate(const double & Bmin, const double & Bmax,
			     Omega_ik * eikonal=NULL);
    virtual double GetValue(const double & B)=0;
    virtual double GetCombinedValue(const double & B)=0;
    virtual double operator()(double B);
    virtual const double & Sigma() const     { return m_sigma; }
    virtual const double & Y() const         { return m_Y; }
    virtual const double & E() const         { return m_E; }
    virtual const double & OriginalY() const { return m_originalY; }
    virtual const double & CutoffY() const   { return m_cutoffY; }
    virtual const double & Accu() const      { return m_accu; }
    virtual const double & XSec() const      { return m_sigma; }
    virtual std::list<Omega_ik *> * Eikonals() { return p_eikonals; }
  };
}
#endif
