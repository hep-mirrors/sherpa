#ifndef METOOLS_Current_Library_VA_0_Novosibirsk4Pi_H
#define METOOLS_Current_Library_VA_0_Novosibirsk4Pi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>

namespace METOOLS { class Total_Width_Base; }

namespace METOOLS {
  // Novosibirsk/CVC 4pi current (Bondar et al., reproduced in
  // tau_two_meson_currents_KS_RChiT.tex Sec."4pi: implementation-
  // ready Novosibirsk/CVC current"). Covers BOTH pi- pi0 pi0 pi0 and
  // pi- pi- pi+ pi0 (a single class, mode determined from flavour
  // content at construction - see FixMode() below).
  //
  // NORMALISATION AND WIDTHS - note a reversal of an earlier decision.
  //
  // This class previously used Sherpa's registered running widths
  // (Total_Width_Base via LineShapes) in place of Bondar's own g_R(s),
  // deliberately and on instruction.  It cannot: Bondar's G tables were
  // fitted against his own width functions, and with registered line
  // shapes the published partial widths come out about 1e-5 too small.
  // The widths therefore come from the parameterisation's own tabulated
  // and analytic functions (Novosibirsk4pi_GTables).
  //
  // POLE MASSES do still come from the particle database, like every
  // other current here.  That costs about 6% against Bondar's own
  // quoted widths - his sigma sits at 0.8 GeV where the database has
  // f0(600) at 0.6 - and is a deliberate trade in favour of one shared
  // registry entry per resonance.
  //
  // The propagators live in Propagator.H (Novo4Pi_Propagator) and the
  // tables and width functions in Novosibirsk4pi_GTables.H; this
  // class is the Lorentz assembly plus the per-mode normalisation.
  //
  class VA_0_Novosibirsk4Pi : public Current_Base {
  private:
    enum class Mode { PiPi0Pi0Pi0, PiPiPiPi0, Unknown };
    Mode m_mode;

    double  m_Mrho, m_Ma1, m_Msigma, m_Momega;
    Total_Width_Base * p_rho_width, * p_a1_width, * p_sigma_width, * p_omega_width;
    Complex m_z;
    double  m_norm;
    // Which piece of the mixed-charge current to return.  Bondar/TAUOLA
    // treat a1pi and omegapi as DISTINCT modes and neglect their
    // interference (note, after Eq.4pi-Jmix-omega): J = J_a1rho +
    // J_a1sigma  OR  J_omegarho.  Adding them coherently, as this class
    // did unconditionally, is a labelled model change.
    //   0 = coherent sum of both (the model change)
    //   1 = a1pi only  (rho + sigma)
    //   2 = omegapi only
    // 1 and 2 exist so the incoherent TAUOLA treatment can be obtained
    // by running them separately and adding the two rates.
    int     m_component;
    // Follow binp.f rather than the paper.  Default ON: the paper alone
    // cannot reproduce the published widths (see the tables header).
    int     m_binp;
    // Bondar's own parameter values (Table 1 of the paper, and the
    // DATA statements of binp.f, which agree).  Deliberately NOT taken
    // from the particle database: this current is a transcription of one
    // specific published-and-coded model, and its G tables were fitted
    // against these values.  Substituting database poles detunes them.
    double  m_ScaleA;   // a1 form-factor scale; this is Lambda^2


    // ---- Bondar's parameterisation ----
    // The propagators live in Propagator.H (Novo4Pi_Propagator) and the
    // width/form-factor functions and tables in
    // Novosibirsk4pi_GTables.H.  This class owns only the Lorentz
    // assembly and the per-mode normalisation below, which is the one
    // piece that is genuinely specific to this current.
    Propagator_Base * p_bRho, * p_bSigma, * p_bA1, * p_bOmega;
    // The per-mode overall factor; no counterpart anywhere in the
    // paper, see the tables header.
    double  Novo4PiModeFactor(const int mode,const double & Q) const;
    ATOOLS::Vec4C t1(const ATOOLS::Vec4D & q1,const ATOOLS::Vec4D & q2,
		     const ATOOLS::Vec4D & q3,const ATOOLS::Vec4D & q4,
		     const ATOOLS::Vec4D & Q) const;
    ATOOLS::Vec4C t2(const ATOOLS::Vec4D & q1,const ATOOLS::Vec4D & q2,
		     const ATOOLS::Vec4D & q3,const ATOOLS::Vec4D & q4,
		     const ATOOLS::Vec4D & Q) const;
    ATOOLS::Vec4C t3(const ATOOLS::Vec4D & q1,const ATOOLS::Vec4D & q2,
		     const ATOOLS::Vec4D & q3,const ATOOLS::Vec4D & q4,
		     const ATOOLS::Vec4D & Q) const;

    void FixMode();

  public:
    VA_0_Novosibirsk4Pi(const ATOOLS::Flavour_Vector& flavs,
			const std::vector<int>& indices,
			const std::string& name);
    ~VA_0_Novosibirsk4Pi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
