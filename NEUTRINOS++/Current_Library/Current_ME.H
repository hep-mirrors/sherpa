#ifndef NEUTRINOS_Current_Library_Current_ME_H
#define NEUTRINOS_Current_Library_Current_ME_H

#include "METOOLS/Main/XYZFuncs.H"
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"

namespace NEUTRINOS {
  class Scatter_Current_Base;

  class Current_ME : public HADRONS::HD_ME_Base {
    Scatter_Current_Base* p_c1, * p_c2;
  public:
    Current_ME(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& decayindices,
               const std::string& name);
    ~Current_ME();
    void SetCurrent1(Scatter_Current_Base* c1) { p_c1 = c1; }
    void SetCurrent2(Scatter_Current_Base* c2) { p_c2 = c2; }
    inline Scatter_Current_Base* GetCurrent1() const { return p_c1; }
    inline Scatter_Current_Base* GetCurrent2() const { return p_c2; }

    double operator()(const ATOOLS::Vec4D_Vector& momenta);
    void   Calculate(const ATOOLS::Vec4D_Vector& p,bool m_anti=false);
    void   Calculate(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F);
  };
}
#endif
