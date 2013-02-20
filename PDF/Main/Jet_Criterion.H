#ifndef PDF_Main_Jet_Criterion_H
#define PDF_Main_Jet_Criterion_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace PDF {

  class Shower_Base;

  struct JetCriterion_Key {
    std::string m_key;
    Shower_Base *p_shower;
    JetCriterion_Key(const std::string &key,
		     Shower_Base *const shower):
      m_key(key), p_shower(shower) {}
  };// end of struct JetCriterion_Key

  class Jet_Criterion {
  public:

    virtual ~Jet_Criterion();

    virtual bool Jets(ATOOLS::Cluster_Amplitude *ampl,int mode=0) = 0;

  };// end of class Jet_Criterion

  typedef ATOOLS::Getter_Function
  <Jet_Criterion,JetCriterion_Key> JetCriterion_Getter;

}// end of namespace PDF

#endif
