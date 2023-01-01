#ifndef BEAM_Main_RelicDensity_Kinematics_H
#define BEAM_Main_RelicDensity_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
  class RelicDensity_Kinematics : public Kinematics_Base {
  private:
    ATOOLS::Info_Key m_sprimekey;
    
    void InitIntegration();
  public:
    RelicDensity_Kinematics(Beam_Base ** beams);
    ~RelicDensity_Kinematics();

    void AssignKeys(ATOOLS::Integration_Info *const info);
    void SetLimits();
    inline const double Sprime() { return m_sprimekey[3]; }
    bool operator()(ATOOLS::Vec4D * moms);
  };
}
#endif
