#ifndef PHASIC_Scales_KFactor_Setter_Base_H
#define PHASIC_Scales_KFactor_Setter_Base_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

namespace PHASIC {

  struct KFactor_Setter_Arguments {
    Process_Base *p_proc;
    std::string m_kfac;
    size_t m_oqcdlo, m_oewlo;
    KFactor_Setter_Arguments(Process_Base *const proc,
			     const std::string &kfac,
			     const size_t &oqcdlo,const size_t &oewlo):
      p_proc(proc), m_kfac(kfac), m_oqcdlo(oqcdlo), m_oewlo(oewlo) {}
  };// end of class KFactor_Setter_Arguments

  class KFactor_Setter_Base {
  public:

    typedef ATOOLS::Getter_Function
    <KFactor_Setter_Base,KFactor_Setter_Arguments> KFactor_Getter_Function;

  protected:

    Process_Base *p_proc;

    bool m_on;

    ATOOLS::Info_Key m_kfkey;

  public:

    KFactor_Setter_Base(Process_Base *const proc);

    virtual ~KFactor_Setter_Base();

    virtual double KFactor() = 0;

    double KFactor2(const bool diff);

    inline Process_Base *Process() const { return p_proc; }

    inline void SetOn(const bool on) { m_on=on; }

    static void ShowSyntax(const size_t i);

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

#endif
