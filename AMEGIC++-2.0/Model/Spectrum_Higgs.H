// modul for parameters

#ifndef Spectrum_Higgs_h
#define Spectrum_Higgs_h

#include "Matrix.H"
#include "Isajet/Isajet.H"
#include "Couplings_EW.H"

namespace AMEGIC {

  class Spectrum_Higgs {
    Couplings_EW CplEW;
    void Masses_LO();
    double v,v1,v2,tanb,mu,M_A0;
    AMATOOLS::Matrix<2> _Z_H;
    AMATOOLS::Matrix<2> _Z_R;
    double M_W,_sin2TW,_ALFA1;
  public:    
    void Interface(Isajet*);
    void Init();
    inline double Z_H(short int i,short int j) {return _Z_H[i][j];}
    inline double Z_R(short int i,short int j) {return _Z_R[i][j];}
    inline double v_1() {return v1;}
    inline double v_2() {return v2;}
    inline double v_0() {return v;}
    inline double tan_b() {return tanb;}
  };
}
#endif
