#ifndef METOOLS_Explicit_CS_Dipole_Kinematics_H
#define METOOLS_Explicit_CS_Dipole_Kinematics_H

#include "METOOLS/Explicit/Dipole_Kinematics.H"


namespace METOOLS {

  class CS_Dipole_Kinematics : public Dipole_Kinematics {
  protected:


  public:

    // constructor
    CS_Dipole_Kinematics(Dipole_Info *const info,Current *const i,
		      Current *const j,Current *const k,
		      Current *const ijt,Current *const kt)
      : Dipole_Kinematics(info,i,j,k,ijt,kt) {}

    // member functions
    void Evaluate();

    void CheckKT2Min();

    void EvaluateFFS(const ATOOLS::Vec4D &p,ATOOLS::Vec4D &q,
                     double &Ai,double &B,double &t);
    void EvaluateVVS(const ATOOLS::Vec4D& p, ATOOLS::Vec4D& q,
                     double &Ai,double &Aj,double &B,double &t);
    void EvaluateFVS(double &Ai,double &t);
    
  };// end of class Dipole_Kinematics
  
}// end of namespace METOOLS

#endif
