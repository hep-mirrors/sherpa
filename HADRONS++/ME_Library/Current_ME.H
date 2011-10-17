#ifndef HADRONS_ME_Library_Current_ME_H
#define HADRONS_ME_Library_Current_ME_H
 
#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class Current_Base;

  class Current_ME : public HD_ME_Base {
    Current_Base* p_c1, * p_c2;
  public:
    Current_ME(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& decayindices,
               const std::string& name);
    ~Current_ME();
    void SetCurrent1(Current_Base* c1);
    void SetCurrent2(Current_Base* c2);
    inline Current_Base* GetCurrent1() const { return p_c1; }
    inline Current_Base* GetCurrent2() const { return p_c2; }
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}
#endif
