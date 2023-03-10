#ifndef HADRONS_ME_Library_Current_ME_H
#define HADRONS_ME_Library_Current_ME_H
 
#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class Current_Base;

  class Current_ME : public HD_ME_Base {
    Current_Base* p_c1, * p_c2;
    double m_factor;
    /** 
     * Select a current from a string using the Getter_Function for HD_Current_Base.
     * 
     * @param current_string String containing the (getter-able) name of the
     * current, and in square brackets the mapping of external particle
     * indices in the decay blob to internal indices in the current.
     * 
     * @return The current it has created from the information. If
     *   it can't find a valid current, it will abort.
     */
    Current_Base* SelectCurrent(ATOOLS::Scoped_Settings s);
    void SetModelParameters(ATOOLS::Scoped_Settings& s);
  public:
    Current_ME(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& decayindices,
               const std::string& name);
    ~Current_ME();
    inline Current_Base* GetCurrent1() const { return p_c1; }
    inline Current_Base* GetCurrent2() const { return p_c2; }
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}
#endif
