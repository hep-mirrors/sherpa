#ifndef Running_AlphaQED_H
#define Running_AlphaQED_H

#include "Flavour.H"
#include "Function_Base.H"

namespace MODEL {
  class Running_AlphaQED : public AMATOOLS::Function_Base {
    const static double m_A[4],m_B[4],m_C[4];
    double m_alpha0,m_alpha_MZ,m_alpha_eff;
    double m_MZ2,m_ecms2,m_split_scale;
  

    double  PiGamma(const APHYTOOLS::Flavour &, double);
  public:
    Running_AlphaQED(const double,const double,const double=0.);

    double operator()(double); 
    double AqedThomson()  { return m_alpha0; }
    double AqedMZ()       { return m_alpha_MZ; }
    double AqedFixed()    { return m_defval; }
    double Aqed(double);
  };

  extern Running_AlphaQED * aqed;

  // according to 
  //   R. Kleiss et al.: CERN yellow report 89-08, vol.3 p.129
  //   Hadronic component from: H. Burkhardt et al.: Z. Phys C43 (89) 497
  // cf. HERWIG / PYTHIA

  //! QED SM Running of Alpha_QED (if swiched on in C-File!)
  /*! cf. "class Couplings_EW" in AMEGIC/Model */
  // fixed-file, fixed-auto, running
  //    std::string    path;

  //    Running_AlphaQED(std::string);
  // initialise alpha0, alpha_MZ and alpha_eff;
  //    void Init(string);

  // running coupling 
  // alpha at Thomson limit (Q=0); used , e.g.,in  splitting functions
  // alpha at MZ
  // alpha at the beam energy (calculated during Init() or read from Parameterfile)
  // possibly running coupling (if switched on)
}
#endif
