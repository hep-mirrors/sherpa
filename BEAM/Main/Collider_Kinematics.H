#ifndef BEAM_Main_Collider_Kinematics_H
#define BEAM_Main_Collider_Kinematics_H

#include "BEAM/Main/Kinematics_Base.H"

namespace BEAM {
  class Collider_Kinematics : public Kinematics_Base {
  private:
    collidermode::code m_mode;

    ATOOLS::Vec4D    m_Pcms, m_fixp_cms[2];
    double           m_S, m_Ecms, m_E1, m_E2;
    double           m_ymin, m_ymax, m_xmin, m_xmax;
    bool             m_asymmetric;
    ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;
    double           m_logtau, m_sprime, m_x[2];

    void CalculateAndSetX(size_t beam,const ATOOLS::Vec4D & p);
    void InitSystem();
    void InitIntegration();
    bool MakeMonochromaticBeams(ATOOLS::Vec4D * moms);
    bool MakeCollinearBeams(ATOOLS::Vec4D * moms);
    void CalculateAndSetX(size_t beam,const ATOOLS::Vec4D * p);
  public:
    Collider_Kinematics(Beam_Base ** beams);
    ~Collider_Kinematics();

    inline const double Sprime() { return m_sprimekey[3]; }
    inline const double Y()      { return m_ykey[2]; }
    void AssignKeys(ATOOLS::Integration_Info *const info);
    void SetLimits();
    bool operator()(ATOOLS::Vec4D * moms);
  };
}
#endif
