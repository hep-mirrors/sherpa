// modul for parameters

#ifndef Spectrum_sNeutrinos_h
#define Spectrum_sNeutrinos_h

#include "Matrix.H"
#include "Isajet/Isajet.H"

namespace AMEGIC {

  class Spectrum_sNeutrinos {
    void Masses_LO();
    double v,tanb,ml[3];
    AMATOOLS::Matrix<3> _Znue;
  public:
    void Init();
    void Interface(Isajet*);
    inline double Znue(short int i,short int j) {return _Znue[i][j];}
  };
}
#endif














