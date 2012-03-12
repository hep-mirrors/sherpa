#ifndef SINIC_EVENT_GENERATOR_BASE_H
#define SINIC_EVENT_GENERATOR_BASE_H
#include "SINIC++/Eikonals/Omega_ik.H"

namespace SINIC {
  class Event_Generator_Base {
  protected:
    Omega_ik * p_eikonal;
    double     m_smin;
  public:
    Event_Generator_Base():p_eikonal(NULL),m_smin(0.){} 
    ~Event_Generator_Base(){};
    virtual Omega_ik * GetEikonal() const {return p_eikonal;}
    virtual const double Smin() const {return m_smin;}
    virtual const bool IsLastRescatter() const { return false; }
  };
}
#endif
