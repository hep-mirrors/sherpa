#ifndef PDF_Main_NLOMC_Base_H
#define PDF_Main_NLOMC_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/Event_Weights.H"
#include "ATOOLS/Phys/NLO_Subevt.H"

#include <string>

namespace ATOOLS {

  class Cluster_Amplitude;
  struct RB_Data;

}

namespace MODEL { class Model_Base; }

namespace PHASIC { class Process_Base; }

namespace PDF {

  class ISR_Handler;
  class Shower_Base;

  class NLOMC_Base {
  protected:

    std::string m_name;

    Shower_Base *p_shower;

    double m_kt2min[2];

    ATOOLS::Event_Weights m_weights;

    int m_subtype;

  public:

    NLOMC_Base(const std::string &name);

    virtual ~NLOMC_Base();

    virtual int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual double KT2(const ATOOLS::NLO_subevt &sub,
		       const double &x,const double &y,
		       const double &Q2) = 0;

    static void ShowSyntax(const int mode);

    inline const std::string &Name() const { return m_name; }

    inline double KT2Min(const int mode) const { return m_kt2min[mode]; }

    inline ATOOLS::Event_Weights Weights() const { return m_weights; }

    inline void SetShower(Shower_Base *const shower) { p_shower=shower; }

    inline int SubtractionType() const { return m_subtype; }

  };// end of class NLOMC_Base

  struct NLOMC_Key {
    MODEL::Model_Base *p_model;
    ISR_Handler *p_isr;
    inline NLOMC_Key(MODEL::Model_Base *const model,
		      ISR_Handler *const isr):
      p_model(model), p_isr(isr) {}
  };//end of struct NLOMC_Key

  typedef ATOOLS::Getter_Function
  <NLOMC_Base,NLOMC_Key> NLOMC_Getter;

}// end of namespace PDF

#endif
