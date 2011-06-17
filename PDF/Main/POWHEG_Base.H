#ifndef PDF_Main_POWHEG_Base_H
#define PDF_Main_POWHEG_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/NLO_Subevt.H"

#include <string>

namespace ATOOLS {

  class Cluster_Amplitude;
  class Data_Reader;
  struct RB_Data;

}

namespace MODEL { class Model_Base; }

namespace PHASIC { class Process_Base; }

namespace PDF {

  class ISR_Handler;
  class Shower_Base;

  typedef std::pair<double,double> ZH_Pair;

  struct Dipole_Value {
    PHASIC::Process_Base *p_proc;
    double m_d, m_b, m_t, m_mu2;
    bool m_trig;
    inline Dipole_Value(PHASIC::Process_Base *const proc=NULL,
			const double &d=0.0,const double &b=0.0,
			const double &t=0.0,const double &mu2=0.0,
			const bool trig=true):
      p_proc(proc), m_d(d), m_b(b), m_t(t), m_mu2(mu2), m_trig(trig) {}
  };// end of struct Dipole_Value

  struct Rho_Map: public std::map<ATOOLS::DDip_ID,Dipole_Value> {
    double m_t0;
    inline Rho_Map(): m_t0(0.0) {}
  };// end of struct Rho_Map

  class POWHEG_Base {
  protected:

    std::string m_name;

    Shower_Base *p_shower;

    double m_kt2min;

  public:

    POWHEG_Base(const std::string &name);

    virtual ~POWHEG_Base();

    virtual Rho_Map GetRho(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual ZH_Pair ZHSplit(const double &B,const double &Qij2,
			    const ATOOLS::RB_Data *rbd) const = 0;

    static void ShowSyntax(const int mode);

    inline const std::string &Name() const { return m_name; }

    inline double KT2Min() const { return m_kt2min; }

    inline void SetShower(Shower_Base *const shower) { p_shower=shower; }

  };// end of class POWHEG_Base

  struct POWHEG_Key {
    MODEL::Model_Base *p_model;
    ISR_Handler *p_isr;
    ATOOLS::Data_Reader *p_read;
    inline POWHEG_Key(MODEL::Model_Base *const model,
		      ISR_Handler *const isr,
		      ATOOLS::Data_Reader *const read):
      p_model(model), p_isr(isr), p_read(read) {}
  };//end of struct POWHEG_Key

  typedef ATOOLS::Getter_Function
  <POWHEG_Base,POWHEG_Key> POWHEG_Getter;

}// end of namespace ATOOLS

#endif
