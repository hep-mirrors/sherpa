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

//! Runcard keys and output labels that distinguish one photon multiplicity from another.
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
  The multiplicity enters only through m_nphotons, which selects a
  Real_Config generated from it, so no n is special.
*/
class Real_Correction : public YFS_Base {
public:

  /*!
    \param pi        process info, already carrying the extra photon legs
    \param nphotons  photons this correction carries (1 = real, 2 = double real, ...)
  */
  Real_Correction(const PHASIC::Process_Info &pi, size_t nphotons);
  virtual ~Real_Correction();

  //! The matrix element at this phase-space point, cuts and comparison applied.
  double Calc_R(const ATOOLS::Vec4D_Vector &p);
  //! The external generator's matrix element, with the symmetry factor applied.
  double Calc_External(const ATOOLS::Vec4D_Vector &p);
  void SoftScan(const ATOOLS::Vec4D_Vector &p);
  //! YFS: IFI false: sum over photon placements of |M|^2 over the coherent |M|^2, from Comix at p
  double NoIFIRatio(const ATOOLS::Vec4D_Vector &p);

  ATOOLS::Cluster_Amplitude *CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const;

  inline bool FailCut() const { return m_failcut; }
  inline void SetProc(PHASIC::Process_Base *proc) { p_proc = proc; p_permproc = NULL; }
  //! one channel (a member of a container group) of the process, with its
  //! flavours in YFS order (the Born channel's legs, then the photons)
  inline void SetChannelProc(PHASIC::Process_Base *proc,
                             const ATOOLS::Flavour_Vector &flavs)
  { p_proc = proc; m_flavs = flavs; p_permproc = NULL; }
  inline size_t NPhotons() const { return m_nphotons; }

  //! Name of the generator that supplies Calc_R's value: the external one
  //! (p_real_me) when set, otherwise the internal one that built p_proc.
  std::string ActiveGenName() const;

  //! Owned: Tree_ME2_Base::GetME2 builds a fresh object for the caller.
  std::unique_ptr<PHASIC::Tree_ME2_Base> p_real_me;   //!< external generator
  std::unique_ptr<PHASIC::Tree_ME2_Base> p_real_me2;  //!< optional second external generator
  MODEL::Coupling_Map m_cpls;

  //! Not owned: the process supplying the internal matrix element, set by SetProc.
  PHASIC::Process_Base *p_proc;
  //! Raw pointer: Cluster_Amplitude is returned to its free list with Delete().
  ATOOLS::Cluster_Amplitude *p_ampl;

  //! <n>_Generator as configured ("Comix"/"Amegic" means no external reference).
  std::string m_gen;

  double m_sym, m_factor;
  bool m_check, m_writemom;
  int m_nmom, m_fill;
  Flavour_Vector m_flavs;

  std::unique_ptr<ME_Compare> p_cmp;

  //! Comix's per-helicity amplitudes for the last point evaluated, cached
  //! when m_keepamps is set; filled in Calc_R (re-entering Differential() elsewhere is unsafe).
  std::vector<METOOLS::Spin_Amplitudes> m_spinamps;
  bool m_keepamps = false;

  /*!
    Comix's per-helicity amplitudes at one arbitrary phase-space point,
    evaluated on demand and cached in m_spinamps. Safe to call from ACRAIC,
    which already runs inside the Born Differential().

    Returns NULL when there is no process, the leg count does not match, or
    Comix produced no amplitudes. The pointer is m_spinamps, valid until the
    next call to this or to Calc_R.
  */
  const std::vector<METOOLS::Spin_Amplitudes> *
  ComixAmplitudes(const ATOOLS::Vec4D_Vector &p,
                  const double propscale = -1.0,
                  const double decscale = -1.0, const size_t deccid = 0);
  const std::vector<METOOLS::Spin_Amplitudes> *
  ComixAmplitudesShifts(const ATOOLS::Vec4D_Vector &p,
                        const PHASIC::Process_Base::Prop_Shifts &shifts);

  //! Nominal weight of the last ComixAmplitudes() evaluation.
  double m_lastcomix = 0.;

  /*!
    YFS order is m_flavs: the Born legs, then the photons. The process's legs
    are in Sherpa's canonical order (a photon before massless quarks).
    m_perm[j] is the process leg of YFS leg j; momenta, masks and spin
    amplitudes go through it.
  */
  ATOOLS::Vec4D_Vector ToProcessOrder(const ATOOLS::Vec4D_Vector &p) const;
  size_t ToProcessMask(size_t mask) const;

protected:

  //! built lazily per process, hence mutable
  mutable std::vector<size_t> m_perm;
  mutable bool m_permid = true;
  //! the process m_perm was built for
  mutable const PHASIC::Process_Base *p_permproc = NULL;
  void UpdatePermutation() const;
  //! the spin amplitudes of the process, relabelled into YFS leg order
  void ToYFSOrder(std::vector<METOOLS::Spin_Amplitudes> &amps) const;

  size_t m_nphotons;
  //! Not owned: points into the static per-multiplicity config table.
  const Real_Config *p_cfg;
  //! Dump momenta/ME values for the external cross-check (Check_Real / Check_RR).
  int m_dump;
  //! Per-instance diagnostic streams.
  std::ofstream m_me_out, m_ps_out;

  static const Real_Config &ConfigFor(size_t nphotons);
};

//! One photon.
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
