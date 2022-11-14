#ifndef BEAM_Main_Collider_Kinematics_H
#define BEAM_Main_Collider_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
class Collider_Kinematics : public Kinematics_Base {
private:
  collidermode m_mode;

  ATOOLS::Vec4D m_p_plus, m_p_minus;
  double m_Ecms;
  double m_ymin, m_ymax, m_xmin, m_xmax;
  ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;
  double m_sprime;

  void InitSystem();
  void InitIntegration() override;
  double CalculateTau() const;

public:
  Collider_Kinematics(std::array<Beam_Base*, 2> beams);
  ~Collider_Kinematics() = default;

  inline double Sprime() override { return m_sprimekey[3]; }
  inline double Y() override { return m_ykey[2]; }
  void AssignKeys(ATOOLS::Integration_Info *const info) override;
  void SetLimits() override;
  bool operator()(ATOOLS::Vec4D_Vector& moms) override;
};
}
#endif
