// modul for parameters

#ifndef Spectrum_sUpquarks_h
#define Spectrum_sUpquarks_h

#include "Matrix.H"
#include "Isajet/Isajet.H"

namespace AMEGIC {

  class Spectrum_sUpquarks {
    void Masses_LO();
    double v,tanb,mu,A_t,mgluino;
    double m[3],mQ[3],mU[3],ws[3],us[3];
    AMATOOLS::Matrix<6> _Zu;
    double _mu;
    AMATOOLS::Matrix<3> _ws,_us;
  public:
    void Init();
    void Interface(Isajet*);
    inline double Zu(short int i,short int j)  {return _Zu[i][j];}
    inline double h()                          {return _mu;}
    inline double w_S(short int i,short int j) {return _ws[i][j];}
    inline double u_S(short int i,short int j) {return _us[i][j];} 
  
 };
}
#endif














