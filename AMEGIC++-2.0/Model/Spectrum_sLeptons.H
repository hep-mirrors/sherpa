// modul for parameters

#ifndef Spectrum_sLeptons_h
#define Spectrum_sLeptons_h

#include "Matrix.H"
#include "Isajet/Isajet.H"

namespace AMEGIC {

  class Spectrum_sLeptons {
    void Masses_LO();
    double v,tanb,mu,A_l;
    double m[3],ml[3],mr[3],ks[3],ls[3];
    AMATOOLS::Matrix<6> _ZLep;
    double _mu;
    AMATOOLS::Matrix<3> _ls,_ks;
    
  public:
    void Init();
    void Interface(Isajet*);
    inline double Z_L(short int i,short int j) {return _ZLep[i][j];}
    inline double h()                          {return _mu;}
    inline double l_S(short int i,short int j) {return _ls[i][j];}
    inline double k_S(short int i,short int j) {return _ks[i][j];}
  };
}
#endif














