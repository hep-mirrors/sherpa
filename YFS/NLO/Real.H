#ifndef YFS_REAL_H
#define YFS_REAL_H


#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/NLO/ME_Compare.H"


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

class Real : public YFS_Base {

public:

  Real(const PHASIC::Process_Info& pi);

  ~Real();
  // void InitME(MODEL::Model_Base *const model);
  double Calc(const ATOOLS::Vec4D_Vector momenta);
   PHASIC::Tree_ME2_Base     *p_real_me;
   PHASIC::Tree_ME2_Base     *p_real_me2;
  double Calc_R(const ATOOLS::Vec4D_Vector& p);
  double Calc_External(const ATOOLS::Vec4D_Vector& p);
  MODEL::Coupling_Map m_cpls;
  PHASIC::Scale_Setter_Base *p_scale;
  PHASIC::Tree_ME2_Base* p_born_me;
  double m_sym, m_factor;
  bool m_check, m_writemom;
  int m_nmom, m_fill;
  Flavour_Vector m_flavs;

  //! Real_Generator as configured (may be "Comix"/"Amegic", meaning "no
  //! external reference, use whatever internal ME_Generator built
  //! p_realproc"). Stored so ActiveGenName() can tell the two apart.
  std::string m_gen;
  //! Name of whichever generator actually supplies Calc_R's value: the
  //! external one (p_real_me, e.g. "OpenLoops") when set, otherwise the
  //! internal ME_Generator that built p_realproc (Comix/Amegic). Use this for
  //! diagnostic labels -- p_realproc->Generator()->Name() alone is wrong
  //! whenever an external Real_Generator is configured, which is the common
  //! case with Compare_Real off.
  std::string ActiveGenName() const;

  //! Explicit phase-space points read from Read_Real_Momenta, if set. Each
  //! entry has m_flavs.size() legs in m_flavs order, incoming legs positive
  //! (i.e. the convention Calc_R's argument `p` and Write_Real_Momenta's
  //! output both use -- CreateAmplitude does its own negation on the way in).
  std::vector<ATOOLS::Vec4D_Vector> m_replaypoints;
  //! Parses a MOMENTA: file (the same format Write_Real_Momenta produces,
  //! and m_check_real's ps-points.yaml dumps already use) into m_replaypoints.
  void ReadMomenta(const std::string &path);
  //! Drives Calc_R once per point in m_replaypoints, in place of normal
  //! phase-space generation, then ends the run.
  void ReplayMomenta();
  //! Appends one phase-space point and its value to the Write_Real_Momenta
  //! files. Called from both of Calc_R's return paths (external-only, and
  //! internal/compared), so the feature works regardless of Compare_Real.
  void WriteMomentum(const ATOOLS::Vec4D_Vector &p, double value);

  ME_Compare *p_cmp;

  PHASIC::Process_Base *p_realproc;
  ATOOLS::Cluster_Amplitude *p_ampl;
  inline bool FailCut() {return m_failcut;}
  inline void SetProc(PHASIC::Process_Base *realproc) {p_realproc=realproc;}
   ATOOLS::Cluster_Amplitude *CreateAmplitude
    (const ATOOLS::Vec4D_Vector &p) const;

  /*!
    Analytic arbiter for the quasi-collinear disagreement between the internal
    and the external real ME. Takes one phase-space point, keeps the hard
    system and the photon energy fixed, and sweeps the photon polar angle
    across the mass angle m/E. In the soft limit the real ME must approach
    S(k) x |M_Born|^2 with the Born independent of the photon direction, so
    ME(theta)/S(theta) must be FLAT - and S is known in closed form,
    mass terms included. Whichever generator's ratio is flat is the correct
    one; no Born, no coupling normalisation and no third generator needed.
  */
  void SoftScan(const ATOOLS::Vec4D_Vector &p);
};
}

#endif