// modul for parameters

#ifndef Spectrum_sDownquarks_h
#define Spectrum_sDownquarks_h

#include "Matrix.H"
#include "Isajet/Isajet.H"

namespace AMEGIC {

  class Spectrum_sDownquarks {
    void Masses_LO();
    double v,tanb,mu,A_b;
    double m[3],mQ[3],mD[3],es[3],ds[3];
    AMATOOLS::Matrix<6> _Zd;
    AMATOOLS::Matrix<3> _es,_ds;
  public:
    void Init();
    void Interface(Isajet*);
    inline double Zd(short int i,short int j) {return _Zd[i][j];}
    inline double e_S(short int i,short int j) {return _es[i][j];}
    inline double d_S(short int i,short int j) {return _ds[i][j];} 

 };
}
#endif














