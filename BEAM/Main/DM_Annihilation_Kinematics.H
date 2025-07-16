#ifndef BEAM_Main_DM_Annihilation_Kinematics_H
#define BEAM_Main_DM_Annihilation_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
  class DM_Annihilation_Kinematics : public Kinematics_Base {
  private:
    ATOOLS::Info_Key m_sprimekey, m_xkey, m_cosxikey;
    double           m_cosxi;

    void InitIntegration();
  public:
    DM_Annihilation_Kinematics(std::array<Beam_Base*, 2> beams);
    ~DM_Annihilation_Kinematics() = default;

    void AssignKeys(ATOOLS::Integration_Info *const info);
    void SetLimits();
    bool operator()(ATOOLS::Vec4D_Vector& p);
  };
}
#endif
