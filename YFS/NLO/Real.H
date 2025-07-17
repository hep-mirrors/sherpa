#ifndef YFS_REAL_H
#define YFS_REAL_H


#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"


namespace PHASIC {
class Tree_ME2_Base;
}

namespace ATOOLS {
class NLO_subevtlist;
}
namespace EXTAMP {
class CS_Dipole;
class Dipole_Wrapper_Process;
}
namespace YFS {
typedef std::map<std::string, ATOOLS::Histogram *> histograms1d;

class Real : public YFS_Base {

public:

  Real(const PHASIC::Process_Info& pi);

  ~Real();
  // void InitME(MODEL::Model_Base *const model);
  double Calc(const ATOOLS::Vec4D_Vector momenta);
   PHASIC::Tree_ME2_Base     *p_real_me;
  double Calc_R(const ATOOLS::Vec4D_Vector& p);
  double Calc_External(const ATOOLS::Vec4D_Vector& p);
  MODEL::Coupling_Map m_cpls;
  PHASIC::Scale_Setter_Base *p_scale;
  PHASIC::Tree_ME2_Base* p_born_me;
  double m_sym, m_factor;
  bool m_check, m_writemom;
  int m_nmom, m_fill;
  std::string m_gen;
  Flavour_Vector m_flavs;
  PHASIC::Process_Base *p_realproc;
  ATOOLS::Cluster_Amplitude *p_ampl;
  histograms1d m_histograms1d;
  inline bool FailCut() {return m_failcut;}
  inline void SetProc(PHASIC::Process_Base *realproc) {p_realproc=realproc;}
   ATOOLS::Cluster_Amplitude *CreateAmplitude
    (const ATOOLS::Vec4D_Vector &p) const;
};
}

#endif