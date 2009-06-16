#ifndef PHASIC_Scales_Tag_Setter_H
#define PHASIC_Scales_Tag_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {

  class Jet_Finder;

  class Tag_Setter: public ATOOLS::Tag_Replacer {
  private:

    Scale_Setter_Base *p_setter;

    ATOOLS::Algebra_Interpreter *p_calc;

    mutable Jet_Finder *p_jf;

    Jet_Finder *JF() const;

  public:
    
    // constructor
    inline Tag_Setter(Scale_Setter_Base *const setter): 
      p_setter(setter), p_calc(NULL), p_jf(NULL) {}
    
    // member functions
    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    // inline functions
    void SetCalculator(ATOOLS::Algebra_Interpreter *const calc) { p_calc=calc; }

  };// end of class Tag_Setter

}// end of namespace PHASIC

#endif
