#ifndef BEAM_Main_Weight_Base_H
#define BEAM_Main_Weight_Base_H

#include "ATOOLS/Math/Vector.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Kinematics_Base.H"
#include <array>

namespace BEAM {
class Weight_Base {
protected:
  Kinematics_Base *p_kinematics;
  std::array<Beam_Base *, 2> p_beams;
  double m_weight;
  std::string m_keyid;

public:
  Weight_Base(Kinematics_Base *kinematics);
  virtual ~Weight_Base();

  virtual void AssignKeys(ATOOLS::Integration_Info *info) = 0;
  virtual bool Calculate(const double &scale) = 0;
  virtual double operator()() {
    return m_weight;
  }
};
} // namespace BEAM
#endif
