
//modul for parameters

#ifndef Spectrum_Charginos_h
#define Spectrum_Charginos_h

#include "Isajet/Isajet.H"
#include "Matrix.H"
#include "Couplings_EW.H"

namespace AMEGIC {

  class Spectrum_Charginos {
    Couplings_EW CplEW;
    void Masses_LO();
    double v,tanb,m_2,mu;
    AMATOOLS::Matrix<2> _Zplus;
    AMATOOLS::Matrix<2> _Zminus;
  public:
    void Init();
    void Interface(Isajet*);
    inline double Zplus(short int i,short int j) {return _Zplus[i][j];}
    inline double Zminus(short int i,short int j) {return _Zminus[i][j];}
  };
}
#endif














