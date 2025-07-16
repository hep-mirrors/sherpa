#ifndef BEAM_Main_RelicDensity_Weight_H
#define BEAM_Main_RelicDensity_Weight_H

#include "BEAM/Main/Weight_Base.H"

namespace BEAM {

  class RelicDensity_Weight : public Weight_Base {
  private:
    ATOOLS::Info_Key m_sprimekey;
    bool        m_relativistic;
    double      m_temperature, m_m[2], m_m2[2], m_w[2];
    long double m_BesselK2[2], m_norm;
  public:
    RelicDensity_Weight(Kinematics_Base * kinematics);
    ~RelicDensity_Weight();

    void AssignKeys(ATOOLS::Integration_Info *const info);
    bool Calculate(const double & scale);
  };
}

#endif
