#ifndef Running_AlphaQED_H
#define Running_AlphaQED_H

//#include <string>
#include "Flavour.H"
#include "Function_Base.H"

namespace APHYTOOLS {
  //! QED SM Running of Alpha_QED (if swiched on in C-File!)
  /*! cf. "class Couplings_EW" in AMEGIC/Model */
  class Running_AlphaQED : public AMATOOLS::Function_Base {
    const static double A[4], B[4], C[4];
    double alpha0,  alpha_MZ,  alpha_eff;
    double splitt_scale;
    int initialised;

    double ecms2;
    // fixed-file, fixed-auto, running
    int mode;          

    //    std::string    path;

    double  PiGamma(const Flavour &, double);
  public:
    //    Running_AlphaQED(std::string);
    Running_AlphaQED();

    // initialise alpha0, alpha_MZ and alpha_eff;
    //    void Init(string);
    void Init();


    double operator()() { return alpha_eff;}

    // running coupling 
    double operator()(double t); 

    // alpha at Thomson limit (Q=0); used , e.g.,in  splitting functions
    double AqedThomson() {return alpha0;}

    // alpha at MZ
    double AqedMZ() {return alpha_MZ;}

    // alpha at the beam energy (calculated during Init() or read from Parameterfile)
    double AqedFixed() {return alpha_eff;}

    // possibly running coupling (if switched on)
    double Aqed(double);
  };

  extern Running_AlphaQED * aqed;
}
#endif
