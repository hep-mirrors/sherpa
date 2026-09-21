#ifndef YFS_REAL_CORRECTION_H
#define YFS_REAL_CORRECTION_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/NLO/ME_Compare.H"
#include <memory>
#include "METOOLS/Main/Spin_Structure.H"

#include <fstream>
#include <string>

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

/*!
  The runcard keys and output labels that distinguish one photon multiplicity
  from another. Everything else about a fixed-order real correction is
  independent of how many photons it carries, which is why Real and RealReal
  were byte-identical copies of each other.
*/
struct Real_Config {
  std::string m_generator;   //!< e.g. "Real_Generator", "RR_Generator"
  std::string m_generator2;  //!< optional second EXTERNAL generator to compare against
  std::string m_compare;     //!< e.g. "Compare_Real", "Compare_RR"
  std::string m_writemom;    //!< e.g. "Write_Real_Momenta"
  std::string m_nmom;        //!< e.g. "N_Real_Momenta"
  std::string m_histdir;     //!< ME_Compare output directory
  std::string m_label;       //!< ME_Compare label
  std::string m_tag;         //!< short tag on the per-point stderr line, e.g. "ROC"
  std::string m_what;        //!< human-readable name for error messages
};

/*!
  One fixed-order real correction of any photon multiplicity.

  This is the single copy of what used to be YFS::Real and YFS::RealReal.
  Those two classes had byte-identical CreateAmplitude() bodies (verified by
  diff after renaming the process pointer) and Calc_External()/Calc_R() bodies
  that differed only in which settings they read and which diagnostics had
  been added to one and not the other. That duplication is how the
  incoming-leg crossing bug came to exist twice and need fixing twice, and how
  the uninitialised external_real read came to exist twice.

  Multiplicity enters ONLY through m_nphotons, which selects a Real_Config -
  and that config is GENERATED from the multiplicity, so no n is special and
  none needs a new class, a new row or a new translation unit.
*/
class Real_Correction : public YFS_Base {
public:

  /*!
    \param pi        process info, already carrying the extra photon legs
    \param nphotons  photons this correction carries (1 = real, 2 = double
                     real, 3 = triple real, ... - no upper bound here; what
                     bounds it is whether an ME provider can be built)
  */
  Real_Correction(const PHASIC::Process_Info &pi, size_t nphotons);
  virtual ~Real_Correction();

  //! The matrix element at this phase-space point, cuts and comparison applied.
  double Calc_R(const ATOOLS::Vec4D_Vector &p);
  //! The external generator's matrix element, with the symmetry factor applied.
  double Calc_External(const ATOOLS::Vec4D_Vector &p);

  ATOOLS::Cluster_Amplitude *CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const;

  /*!
    Analytic arbiter for the quasi-collinear disagreement between the internal
    and the external real ME. Takes one phase-space point, keeps the hard
    system and the photon energy fixed, and sweeps the photon polar angle
    across the mass angle m/E. In the soft limit the real ME must approach
    S(k) x |M_Born|^2 with the Born independent of the photon direction, so
    ME(theta)/S(theta) must be FLAT - and S is known in closed form,
    mass terms included. Whichever generator's ratio is flat is the correct
    one; no Born, no coupling normalisation and no third generator needed.

    Single-photon only so far: the sweep is over ONE photon's polar angle.
    Generalising it is the natural replacement for the cross-generator check
    at n=3, where an external double-real-plus-one is unlikely to exist.
  */
  void SoftScan(const ATOOLS::Vec4D_Vector &p);

  inline bool FailCut() const { return m_failcut; }
  inline void SetProc(PHASIC::Process_Base *proc) { p_proc = proc; }
  inline size_t NPhotons() const { return m_nphotons; }

  //! Name of whichever generator actually supplies Calc_R's value: the
  //! external one (p_real_me, e.g. "OpenLoops") when set, otherwise the
  //! internal ME_Generator that built p_proc (Comix/Amegic). Use this for
  //! diagnostic labels -- p_proc->Generator()->Name() alone is wrong whenever
  //! an external <n>_Generator is configured, which is the common case with
  //! Compare_Real/Compare_RR/... off.
  std::string ActiveGenName() const;

  //! Owned: Tree_ME2_Base::GetME2 builds a fresh object for the caller.
  std::unique_ptr<PHASIC::Tree_ME2_Base> p_real_me;   //!< external generator
  std::unique_ptr<PHASIC::Tree_ME2_Base> p_real_me2;  //!< optional second external generator
  MODEL::Coupling_Map m_cpls;

  //! Not owned: the Sherpa process supplying the internal (Comix/Amegic)
  //! matrix element, set by SetProc.
  PHASIC::Process_Base *p_proc;
  /*!
    Deliberately a raw pointer: Cluster_Amplitude comes from a free list and is
    returned to it with its own Delete(), never with operator delete, so no
    standard smart pointer owns it correctly.
  */
  ATOOLS::Cluster_Amplitude *p_ampl;

  //! <n>_Generator as configured (may be "Comix"/"Amegic", meaning "no
  //! external reference, use whatever internal ME_Generator built p_proc").
  //! Stored so ActiveGenName() can tell the two apart.
  std::string m_gen;

  double m_sym, m_factor;
  bool m_check, m_writemom;
  int m_nmom, m_fill;
  Flavour_Vector m_flavs;

  std::unique_ptr<ME_Compare> p_cmp;

  /*!
    Comix's per-helicity amplitudes for the LAST point evaluated, cached when
    m_keepamps is set. Filled in Calc_R right after the Differential() that
    computes them, which is the only safe place: re-entering Differential()
    from elsewhere in the event (e.g. from CEEX) segfaults.

    Measured to reproduce the matrix element exactly -
    Differential/SumSquare = 1/(4 n!) to the last bit over 5864 points - so
    these are a drop-in source for the CEEX spin structure.
  */
  std::vector<METOOLS::Spin_Amplitudes> m_spinamps;
  bool m_keepamps = false;
  inline void KeepAmplitudes(bool k) { m_keepamps = k; }

protected:

  size_t m_nphotons;
  //! Not owned: points into the static per-multiplicity config table.
  const Real_Config *p_cfg;
  //! Whether to dump momenta/ME values for the external cross-check, from
  //! Check_Real / Check_RR according to multiplicity.
  int m_dump;
  //! Per-instance diagnostic streams. These were file-scope globals in each
  //! of the two former translation units, which is safe only while there is
  //! exactly one instance of each multiplicity.
  std::ofstream m_me_out, m_ps_out;

  static const Real_Config &ConfigFor(size_t nphotons);
};

//! One photon. Thin alias so existing call sites and member names are unchanged.
class Real : public Real_Correction {
public:
  Real(const PHASIC::Process_Info &pi) : Real_Correction(pi, 1) {}
};

//! Two photons.
class RealReal : public Real_Correction {
public:
  RealReal(const PHASIC::Process_Info &pi) : Real_Correction(pi, 2) {}
};

}

#endif
