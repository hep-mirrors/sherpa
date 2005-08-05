#ifndef Sudakov_h
#define Sudakov_h

#include "Splitting_Function_Group.H"
#include "Parton.H"
#include "Running_AlphaS.H"

namespace CS_SHOWER {
  class Sudakov : public Splitting_Function_Group {
    double m_k0sq, m_asmax, m_nspect;
    double m_deltaz, m_kperp2, m_z, m_y, m_scale;
    MODEL::Running_AlphaS * p_as;

    void ProduceT();
    bool Veto(double,double);
    bool KinCheck(double,double);
    bool Splitting();
    bool Coupling();
  public:
    Sudakov();
    ~Sudakov();

    bool   Dice(Parton *,Parton *,const double=0.);
    double OverIntegrated(const double,const double,const double);
  };
}

#endif
