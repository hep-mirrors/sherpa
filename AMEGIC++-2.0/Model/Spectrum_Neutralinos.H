// modul for parameters

#ifndef Spectrum_Neutralinos_h
#define Spectrum_Neutralinos_h

#include "Matrix.H"
#include "Isajet/Isajet.H"
#include "Couplings_EW.H"

namespace AMEGIC {

  class Spectrum_Neutralinos {
    Couplings_EW CplEW;
    void Masses_LO();
    double v,tanb,m_2,m_3,mu;
    AMATOOLS::Matrix<4> _Z_N;
  public:
    void Init();
    void Interface(Isajet*);
    inline double Z_N(short int i,short int j) {return _Z_N[i][j];}
  };
}
#endif














