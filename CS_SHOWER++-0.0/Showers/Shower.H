#ifndef Shower_h
#define Shower_h

#include "ISR_Handler.H"
#include "Sudakov.H"
#include "Singlet.H"
#include "Kinematics_Base.H"

namespace CS_SHOWER {
  class Shower {
  private:
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    ATOOLS::Flavour    m_flavB, m_flavC;
    Singlet      *     p_actual;
    All_Singlets *     p_all;
    Sudakov            m_sudakov;
    
    PLiter SelectSplitting(double &);
    bool   EvolveSinglet(Singlet *);
    bool   PerformSplitting();
 public:
    Shower(PDF::ISR_Handler * isr);
    ~Shower();
    bool   EvolveShower(All_Singlets *);
  };
}

#endif
