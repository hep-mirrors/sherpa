#ifndef SHERPA_SoftPhysics_Hadron_Rescattering_Handler_H
#define SHERPA_SoftPhysics_Hadron_Rescattering_Handler_H

#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"

namespace SHERPA {
  class Hadron_Rescattering_Handler {
  private:
    HADRON_RESCATTERING::Hadron_Rescatterings  * p_rescattering;
    std::string m_name;
  public :
    Hadron_Rescattering_Handler();
    ~Hadron_Rescattering_Handler();
    void CleanUp(const size_t & mode=0) {}

    void Output() {}
    inline const std::string & Name() const { return m_name; }
  };
}

#endif

