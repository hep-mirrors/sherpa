#ifndef PHASIC_Scales_Tag_Setter_H
#define PHASIC_Scales_Tag_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {

  class Tag_Setter: public ATOOLS::Tag_Replacer {
  private:

    Scale_Setter_Base *p_setter;

    ATOOLS::Algebra_Interpreter *p_calc;

  public:
    
    // constructor
    inline Tag_Setter(Scale_Setter_Base *const setter): 
      p_setter(setter), p_calc(NULL) {}
    
    // member functions
    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

    void SetTags(ATOOLS::Algebra_Interpreter *const calc);

    // inline functions
    inline void SetCalculator(ATOOLS::Algebra_Interpreter *const calc)
    { p_calc=calc; }

  };// end of class Tag_Setter

}// end of namespace PHASIC

#endif
