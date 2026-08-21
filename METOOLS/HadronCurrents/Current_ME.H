#ifndef METOOLS_HadronCurrents_Current_ME_H
#define METOOLS_HadronCurrents_Current_ME_H

#include "METOOLS/HadronCurrents/ME_Base.H"
#include "METOOLS/HadronCurrents/Current_Base.H"

  
namespace METOOLS {
  class Current_ME : public ME_Base {
    double m_factor;
    Current_Base * p_c1, * p_c2;
    Current_Base * SelectCurrent(ATOOLS::Scoped_Settings s);
    void SetModelParameters(ATOOLS::Scoped_Settings& s);
  public:
    Current_ME(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& decayindices,
               const std::string& name);
    ~Current_ME();
    inline Current_Base * GetCurrent1() const { return p_c1; }
    inline Current_Base * GetCurrent2() const { return p_c2; }
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}
#endif
