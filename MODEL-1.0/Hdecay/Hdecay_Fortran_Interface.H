#ifndef Hdecay_Fortran_Interface_h
#define Hdecay_Fortran_Interface_h

#include "Spectrum_Generator_Base.H"
#include "Data_Read.H"

namespace HDECAY {
  class Hdecay_Fortran_Interface : public MODEL::Spectrum_Generator_Base {
  private :
    double hmass;
    double brb,brtau,brmu,brc,brs,brt;
    double brw,brz,brgg,brgamma,brzgamma,hwidth;

    void   CalculateEffectiveCouplings(std::string);
  protected :

  public :
    Hdecay_Fortran_Interface(AORGTOOLS::Data_Read *,MODEL::Model_Base *);
    ~Hdecay_Fortran_Interface() {};
    void Run(std::string);
    void FillMasses() {}
    void FillDecays() {}
  };
}

#endif


