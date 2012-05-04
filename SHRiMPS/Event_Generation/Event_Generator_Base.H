#ifndef SHRIMPS_EVENT_GENERATOR_BASE_H
#define SHRIMPS_EVENT_GENERATOR_BASE_H
#include "SHRiMPS/Eikonals/Omega_ik.H"

namespace SHRIMPS {
  class Event_Generator_Base {
  protected:
    Omega_ik * p_eikonal;
    double     m_smin;
  public:
    Event_Generator_Base(): p_eikonal(NULL),m_smin(0.) {} 
    ~Event_Generator_Base() {};
    virtual Omega_ik * GetEikonal() const { return p_eikonal;}
    virtual double Smin()           const { return m_smin;}
    virtual bool IsLastRescatter()  const { return false; }
  };
}
#endif
