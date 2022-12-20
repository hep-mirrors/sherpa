#ifndef BEAM_Main_Collider_Kinematics_H
#define BEAM_Main_Collider_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
class Collider_Kinematics : public Kinematics_Base {
private:
  collidermode::code m_mode;

  ATOOLS::Vec4D m_p_plus, m_p_minus;
  double m_Ecms;
  double m_ymin, m_ymax, m_xmin, m_xmax;
  ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;
  double m_sprime;

  void InitSystem();
  void InitIntegration();
  double CalculateTau();

public:
  Collider_Kinematics(Beam_Base **beams);
  ~Collider_Kinematics();

  inline const double Sprime() { return m_sprimekey[3]; }
  inline const double Y() { return m_ykey[2]; }
  void AssignKeys(ATOOLS::Integration_Info *const info);
  void SetLimits();
  bool operator()(ATOOLS::Vec4D_Vector& moms);
};
}
#endif
