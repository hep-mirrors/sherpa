#ifndef Isajet_h
#define Isajet_h

#include "Matrix.H"

namespace AMEGIC { 
  class Isajet {
  public:
    Isajet();
    void Chargino(AMATOOLS::Matrix<2>&,AMATOOLS::Matrix<2>&);
    void Neutralino(AMATOOLS::Matrix<4>&);
    void sNeutrino(AMATOOLS::Matrix<3>&);
    void Higgs(AMATOOLS::Matrix<2>&,AMATOOLS::Matrix<2>&,double&);
    void sUps(AMATOOLS::Matrix<6>&,double&,AMATOOLS::Matrix<3>&,AMATOOLS::Matrix<3>&);
    void sDowns(AMATOOLS::Matrix<6>&,AMATOOLS::Matrix<3>&,AMATOOLS::Matrix<3>&);
    void sLeptons(AMATOOLS::Matrix<6>&,double&,AMATOOLS::Matrix<3>&,AMATOOLS::Matrix<3>&);
    void Gluino(double&);
  };
}
#endif
