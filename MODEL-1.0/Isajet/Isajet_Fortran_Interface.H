#ifndef Isajet_Fortran_Interface_h
#define Isajet_Fortran_Interface_h

#include "Spectrum_Generator_Base.H"
#include "Data_Read.H"
#include "Matrix.H"

namespace ISAJET { 
  class Isajet_Fortran_Interface : public MODEL::Spectrum_Generator_Base {
    void CharginoMasses();
    void NeutralinoMasses();
    void GluinoMasses();
    void HiggsMasses();
    void sNeutrinoMasses();
    void sUpMasses();
    void sDownMasses();
    void sLeptonMasses();

  public:
    Isajet_Fortran_Interface(AORGTOOLS::Data_Read *,MODEL::Model_Base *);
    ~Isajet_Fortran_Interface();
    void Run(std::string);
    void FillMasses();
    void FillDecays();
  };
}
#endif
