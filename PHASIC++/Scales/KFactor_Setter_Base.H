#ifndef PHASIC_Scales_KFactor_Setter_Base_H
#define PHASIC_Scales_KFactor_Setter_Base_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

namespace PHASIC {

  struct KFactor_Setter_Arguments {
    Process_Base *p_proc;
    std::string m_kfac;
    KFactor_Setter_Arguments(const std::string &kfac):
      m_kfac(kfac) {}
  };// end of class KFactor_Setter_Arguments

  class KFactor_Setter_Base {
  public:

    typedef ATOOLS::Getter_Function
    <KFactor_Setter_Base,KFactor_Setter_Arguments> KFactor_Getter_Function;

  protected:

    Process_Base *p_proc;

    double m_weight;

    bool m_on;

  public:

    KFactor_Setter_Base(const KFactor_Setter_Arguments &args);

    virtual ~KFactor_Setter_Base();

    virtual double KFactor() = 0;

    inline Process_Base *Process() const { return p_proc; }

    inline void SetOn(const bool on) { m_on=on; }

    static void ShowSyntax(const size_t i);

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

#endif
