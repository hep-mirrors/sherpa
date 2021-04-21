#ifndef PHASIC_Selectors_Fastjet_Selector_Base_H
#define PHASIC_Selectors_Fastjet_Selector_Base_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"

namespace PHASIC {

  class Fastjet_Selector_Base : public Selector_Base {

  public:
    Fastjet_Selector_Base(const std::string& name,
                          Process_Base* const proc,
                          ATOOLS::Scoped_Settings s);
    ~Fastjet_Selector_Base();

  protected:
    double m_ptmin, m_etmin, m_delta_r, m_f, m_eta, m_y;
    int m_nj, m_eekt;
    fjcore::JetDefinition * p_jdef;

  public:
    static void PrintCommonInfoLines(std::ostream& str, size_t width);
  };

}

#endif
