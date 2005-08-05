#ifndef Shower_h
#define Shower_h

#include "Sudakov.H"
#include "Singlet.H"
#include "Kinematics_Base.H"

namespace CS_SHOWER {
  class Shower {
  private:
    Sudakov         m_sudakov;
    Kinematics_FF   m_kin;
    ATOOLS::Flavour m_flavB, m_flavC;
    Singlet      *  p_actual;
    All_Singlets *  p_all;

    PLiter SelectSplitting(double &);
    bool   EvolveSinglet(Singlet *);
    bool   PerformSplitting();
 public:
    Shower();
    ~Shower();
    
    bool   EvolveShower(All_Singlets *);
  };
}

#endif
