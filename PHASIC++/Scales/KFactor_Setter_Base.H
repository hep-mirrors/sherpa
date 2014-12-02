#ifndef PHASIC_Scales_KFactor_Setter_Base_H
#define PHASIC_Scales_KFactor_Setter_Base_H

#include "ATOOLS/Org/Exception.H"
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

    // compute KFactor according to kinematics and flavours of p_proc
    virtual double KFactor() = 0;
    // compute KFactor according to kinematics and flavours of NLO subevent
    virtual double KFactor(const ATOOLS::NLO_subevt& evt) {
      THROW(not_implemented, "No NLO (R-S) version of this method implemented");
    }

    inline Process_Base *Process() const { return p_proc; }

    inline void SetOn(const bool on) { m_on=on; }

    static void ShowSyntax(const size_t i);

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

#endif
