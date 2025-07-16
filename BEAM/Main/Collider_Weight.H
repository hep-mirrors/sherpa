#ifndef BEAM_Main_Collider_Weight_H
#define BEAM_Main_Collider_Weight_H

#include "BEAM/Main/Weight_Base.H"

namespace BEAM {
class Collider_Weight : public Weight_Base {
private:
  collidermode m_mode;
  ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;

public:
  Collider_Weight(Kinematics_Base *kinematics);
  ~Collider_Weight();

  void AssignKeys(ATOOLS::Integration_Info *const info) override;
  bool Calculate(const double &scale) override;
  double operator()() override;
};
} // namespace BEAM
#endif
