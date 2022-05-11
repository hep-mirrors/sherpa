#ifndef PHASIC_Scales_Variable_Scale_Setter_H
#define PHASIC_Scales_Variable_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/Tag_Setter.H"

namespace PHASIC {

  class Core_Scale_Setter;

  class Variable_Scale_Setter: public Scale_Setter_Base {
  protected:

    Core_Scale_Setter *p_core;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

  public:

    Variable_Scale_Setter(const Scale_Setter_Arguments &args);

    ~Variable_Scale_Setter();

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode) override;

    PDF::CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const override;

  private:
    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);


  };// end of class Scale_Setter_Base

}// end of namespace PHASIC


#endif
