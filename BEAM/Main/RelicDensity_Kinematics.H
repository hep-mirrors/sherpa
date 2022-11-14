#ifndef BEAM_Main_RelicDensity_Kinematics_H
#define BEAM_Main_RelicDensity_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
  class RelicDensity_Kinematics : public Kinematics_Base {
  private:
    ATOOLS::Info_Key m_sprimekey;

    void InitIntegration();
  public:
    RelicDensity_Kinematics(std::array<Beam_Base*, 2> beams);
    ~RelicDensity_Kinematics() = default;

    void AssignKeys(ATOOLS::Integration_Info *const info);
    void SetLimits();
    inline double Sprime() { return m_sprimekey[3]; }
    bool operator()(ATOOLS::Vec4D_Vector& moms);
  };
}
#endif
