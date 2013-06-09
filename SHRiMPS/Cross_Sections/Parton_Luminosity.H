#ifndef SHRIMPS_Cross_Sections_Parton_Luminosity_H
#define SHRIMPS_Cross_Sections_Parton_Luminosity_H

#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Parton_Luminosity {

    class Kernel : public ATOOLS::Function_Base {
    private:
      double          m_smin,m_smax,m_eta,m_shat,m_maxdl;
      Continued_PDF * p_pdf[2];
    public:
      Kernel() {}
      Kernel(const double & smin,const double & smax) : 
	m_smin(smin), m_smax(smax), m_eta(0.), m_shat(0.) {}
      inline ~Kernel() {}
      double operator()(double y);
      inline double Weight(const double & shat) const {
	return pow(shat/m_smin,m_eta); 
      }
      inline double Weight(const double & shat,const double & eta) const {
	return pow(shat/m_smin,eta); 
      }
      inline void SetISR(Continued_PDF *pdf[2]) {
	p_pdf[0] = pdf[0];
	p_pdf[1] = pdf[1];
      }
      inline void Reset()                         { m_maxdl = 0.; }
      inline void SetExponent(const double & eta) { m_eta   = eta;  }
      inline void SetShat(const double & shat)    { m_shat  = shat; }
      inline void SetSmin(const double & smin)    { m_smin  = smin; }
      inline const double & MaxDL() const         { return m_maxdl; }
    };

  private:
    bool   m_adjust;
    double m_Emin,m_Ecms,m_smin,m_smintest,m_smax;
    int    m_shatsteps;
    double m_deltalog,m_maxdl;

    std::map<Omega_ik *,double, eikcomp> m_deltalogs, m_smins, m_maxdls;
    Kernel                       m_kernel;
    Continued_PDF              * p_pdf[2];

    Omega_ik * p_eikonal;
    double     m_b1,m_b2,m_eta;

    int    m_test;
    void   Tests(Omega_ik * eikonal);
    double CalculateTotalXSec(double & maxdl);
  public:
    Parton_Luminosity() {}
    Parton_Luminosity(const double & Emin,const double & Ecms,
		      const int & shatsteps=1000,const int & test=0);
    ~Parton_Luminosity();
    void   FillGrids(std::list<Omega_ik *> * eikonals);
    void   SetEikonal(Omega_ik * eikonal);
    inline double Weight(const double & s) const{
      return m_kernel.Weight(s);
    }
    inline double Weight(const double & s,const double & eta) const{
      return m_kernel.Weight(s,eta);
    }
    inline void SetPDF(Continued_PDF *const pdf,const int & beam) {
      if (beam!=0 && beam!=1) return;
      p_pdf[beam] = pdf;
    }
    inline const double & Smin() const { return m_smin; }
  };
}

#endif
