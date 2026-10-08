#ifndef METOOLS_Current_Library_FormFactor_Base_H
#define METOOLS_Current_Library_FormFactor_Base_H

#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>

namespace METOOLS {
  enum class ff_model {
    none      = 0,
    point     = 1,
    dipole    = 2,
    // --- 100s: Kuehn-Santamaria family (coherent running-width
    //     Breit-Wigner sums) - consecutive numbering, one entry per
    //     propagator/lineshape variant. ---
    KS        = 100,
    // Gounaris-Sakurai vector-meson propagator variant of KS, cf.
    // arXiv:1509.09140 Eq.(2.2), used e.g. by Belle/ALEPH/CLEO. See
    // Propagator.[C,H] (resonance_type::GS) for the actual math. This
    // is a KS-family model (same resonance content/weights, different
    // propagator functional form), hence grouped in the 100s rather
    // than kept as an unrelated standalone number.
    GS        = 101,
    // Second, distinct "CLEO/TAUOLA-type" Kpi/Keta resonance content:
    // K*(892)+w*K*(1680) (vs. KS=100's FM95 K*(892)+beta*K*(1410)) -
    // previously toggled via a hidden "Kpi_CLEO_tune"/"Keta_CLEO_tune"
    // model parameter under the SAME model number (100); now its own
    // selectable family member, per request. Only affects Kpi_plus/
    // Keta_plus for now - requesting this for pipi_plus/KK_plus (which
    // have no such second tune implemented) falls back to KS=100 (see
    // the family-fallback policy below).
    KS_CLEO   = 102,
    // KS_CLEO's Kpipi variant, plus an added f0(500)/sigma admixture
    // in the T_rho^(1) (pipi sub-system) tower - neither FM95 nor the
    // CLEO K1 alternative include any scalar contribution there
    // originally; this is a NEW, requested addition (not from a
    // specific literature source) reusing Sherpa's own registered
    // f0(600) running width. See F1_0_KPiPi::Construct().
    KS_f0     = 103,
    CLEO_f0   = 104,
    // (105 PREVIOUSLY held a Kiers-et-al. CLEO+F3 model, arXiv:0808.1707.
    //  Removed: measured WORSE than plain CLEO on BaBar 2007 and its F3
    //  normalisation was never calibrated - their N/N3 are in a convention
    //  with no explicit CKM factor, so their numbers cannot be used directly.
    //  See the note in FS_0_KPiPi::FF_KS.  It was then briefly "GS_Belle", the
    //  Belle 2008 pion form-factor fit with its numbers hardcoded - which was
    //  just the Gounaris-Sakurai form (101) carrying one specific fit.  Since
    //  model numbers mean the FUNCTIONAL FORM and local ParameterSets carry
    //  the fit, that is FORM_FACTOR = 101 with Resonance_Parameters = local
    //  and ParameterSets/Belle2008_Fujikawa.yaml.  RETIRED 2026-09-22: a
    //  config naming 105 now fails in FormFactor_Base's constructor with that
    //  instruction rather than silently selecting something else.)
    retired_105 = 105,
    KS_flatte = 106,
    // 107 is a TUNE OF OURS, not a published model, and deliberately has no
    // ParameterSets file: it is the CLEO/Mileo K pi pi parameterisation
    // (Mileo, Kiers & Szynkman, arXiv:1410.1909, itself completing CLEO's
    // Asner et al. PRD 62 072006) with the constant-width K1(1270) replaced
    // by the two-channel Flatte of model 106.
    //
    // Why: the K1(1270) is not one resonance.  In the chiral unitary approach
    // it has TWO poles - roughly 1195/246 MeV coupling to K* pi, and 1284/146
    // coupling to K rho - and those are exactly the two sub-chains of this
    // current.  A single Breit-Wigner cannot represent that, which is the
    // likeliest reason every analysis fits a different K1(1270) width (FM95
    // 0.090, CLEO 0.26, the registry 0.20).  Kuhn, Mirkes and Willibald
    // (hep-ph/9712263) traced their own K pi pi disagreement to "the K1
    // widths".  A Flatte is the two-channel amplitude, so it is a cruder
    // version of the right object.
    //
    // CONSEQUENCE, and why this cannot be a parameter set: A, C and N were
    // fixed against the constant-width K1(1270).  Changing that line shape
    // invalidates them - A and C are that resonance's coefficients, and N
    // normalises a width computed with it.  They are SEEDED from Mileo's
    // values here and MUST be re-determined against data before any rate from
    // 107 is quoted.  Same situation as xi in 106, which had to move from
    // FM95's 0.33 to 0.70 once the K1(1270) became a Flatte.
    CLEO_flatte = 107,
    // --- 200s: Resonance Chiral Theory / chiral-resonance family
    //     (tree-level exchange, full RChL currents, and dispersive
    //     refinements built on RChL input) - consecutive numbering. ---
    RChiPT    = 200,
    // Full Resonance Chiral Lagrangian currents of arXiv:1203.3955
    // (base two-/three-meson currents) with the sigma-meson extension
    // of arXiv:1310.1053 for the 3pi channel. See RChL_Functions.[C,H]
    // for the shared building blocks, and FF_0_PP.C / FF_0_PPP.C for
    // the per-channel assembly. Select via FORM_FACTOR = 201.
    RChL2012  = 201,
    // Combined dispersive(low-s)/modified-RChL(high-s) parametrization
    // of arXiv:1509.09140 Eq.(2.3)ff (Gomez Dumm & Roig). NOT YET
    // IMPLEMENTED: the low-s branch needs a numerical dispersion
    // integral over the pi pi P-wave phase shift, and the high-s
    // branch needs the two-meson loop functions B22^pi, B22^K, neither
    // of which exist in this code yet. Selecting this currently throws
    // (see FF_0_PP_Base::operator()). Kept as an enum value so callers/
    // config files can already reference it and get a clear error
    // rather than silently falling through to ff_model::unknown.
    combRChL  = 202,
    // --- 300s: dispersive family -- the two-body dynamics come from a
    //     solved dispersion relation rather than from line shapes. ---
    // Khuri-Treiman representation of gamma* -> 3 pi with two
    // subtractions, tabulated in KT3Pi_Tables.[C,H].  Implemented by
    // FF_0_Isoscalar3Pi only, and only for sqrt(q^2) in 0.740-1.060 GeV;
    // see that class for why the window.
    KT        = 300,
    unknown   = 999
  };

  // Family-fallback policy (request #2, see FF_0_PP_Base/FF_0_PPP_Base
  // ::operator()): if the SPECIFICALLY requested model has no
  // dedicated construction for a given channel, fall back to the
  // "family base" model (KS=100 for any 100s model, RChiPT=200 for
  // any 200s model) if THAT is implemented for the channel, and only
  // fall back further to the constant (isospin/CKM-only) form factor
  // if even the family base isn't available either. KS=100 and
  // RChiPT=200 themselves have no further family fallback (they ARE
  // the base).
  inline bool IsKSFamily(const ff_model & m)     { return int(m)>=100 && int(m)<200; }
  inline bool IsRChiPTFamily(const ff_model & m) { return int(m)>=200 && int(m)<300; }

  class GeneralModel;
  
  struct FF_Parameters {
    ff_model                       m_ffmodel;
    const ATOOLS::Flavour_Vector & m_flavs;
    std::vector<double>            m_masses, m_masses2;
    std::vector<int>               m_pi;
    std::map<std::string,double>   m_parameterMap;
    std::string                    m_name;
    GeneralModel                 * p_model;

    FF_Parameters(const ff_model & ffmodel,
		  const ATOOLS::Flavour_Vector & flavs,
		  const std::vector<int> & decayindices,
		  std::map<std::string,double> & parameterMap,
		  const std::string & name,
		  GeneralModel * model=NULL) :
      m_ffmodel(ffmodel), m_flavs(flavs), m_pi(decayindices),
      m_parameterMap(parameterMap), m_name(name),
      p_model(model)
    {
      m_masses.resize(m_flavs.size());
      m_masses2.resize(m_flavs.size());      
      for (size_t i=0;i<m_flavs.size();i++) {
	m_masses[i]  = m_flavs[i].Mass(true);
	m_masses2[i] = m_masses[i]*m_masses[i]; 
      }
    }
  };
  
  class FormFactor_Base {
  protected:
    ff_model               m_ffmodel;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<double>    m_masses, m_masses2;
    std::vector<int>       m_pi;
    std::string            m_name;
    GeneralModel         * p_model;
  public:
    FormFactor_Base(const FF_Parameters & params);
    inline virtual ~FormFactor_Base() {};
    virtual Complex operator()(const ATOOLS::Vec4D_Vector& moms) = 0;
    inline void        SetName(std::string name) { m_name = name; }
    inline std::string Name() const              { return m_name; }
  };
}

#define DECLARE_FF_GETTER(CLASS,TAG)				\
  DECLARE_GETTER(CLASS,TAG,METOOLS::FormFactor_Base,METOOLS::FF_Parameters); \
  void ATOOLS::Getter<METOOLS::FormFactor_Base,METOOLS::FF_Parameters,CLASS>:: \
  PrintInfo(std::ostream &str,const size_t width) const { str<<"FF::"<<TAG; }

typedef ATOOLS::Getter_Function
<METOOLS::FormFactor_Base,METOOLS::FF_Parameters> FF_Getter;

#endif
