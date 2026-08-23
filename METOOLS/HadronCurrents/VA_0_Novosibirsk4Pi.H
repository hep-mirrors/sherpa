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
  // IMPORTANT DEPARTURE FROM THE SOURCE, per explicit instruction:
  // wherever Sherpa already has a genuine, physically-motivated
  // registered running width (Total_Width_Base, via LineShapes),
  // THAT is used in place of Bondar's own analytic g_R(s) functions -
  // this affects all four propagators D_rho, D_a1, D_sigma, D_omega.
  // This is a deliberate, non-trivial choice:
  //   - D_rho, D_omega: straightforward substitutions (Sherpa's
  //     rho(770)+/omega(782) lineshapes are themselves built from
  //     genuine rho->pipi / omega->rhopi->3pi decay channels).
  //   - D_sigma: uses Sherpa's F_0_500_Lineshape (kf_f_0_600).
  //   - D_a1: this ELIMINATES the single biggest implementation
  //     blocker in the source note - Bondar's own g_a1(s) is DEFINED
  //     as a 3-body phase-space integral of the a1->rho pi + sigma pi
  //     decay amplitude itself (a genuinely self-referential
  //     calculation the note says is normally handled by pre-
  //     tabulating and interpolating). Sherpa's OWN registered
  //     A1_1260_plus_Lineshape ALREADY implements essentially this
  //     same kind of calculation (a1->{rho(770)+rho(1450)}pi and
  //     a1->{f0(600)+f0(1370)}pi, via the V_PPP/V_PPP_scalar 3-body
  //     partial-width machinery, pre-integrated over a 500-point
  //     grid) - reusing it directly avoids re-deriving this from
  //     scratch, at the cost of no longer numerically matching
  //     Bondar's own specific g_a1(s) formula exactly.
  // Pole masses M_R are ALSO taken from Sherpa's own registered
  // Flavour().HadMass() for each resonance (not Bondar's own specific
  // table values, e.g. Mrho=0.7761 GeV) - same "reuse everywhere
  // possible" consistency choice, small numerical differences
  // expected as a result.
  //
  // The a1 form factor F_a1(q) (Eq.4pi-Fa1, scale Lambda) and the
  // sigma-pi/rho-pi amplitude ratio z (Eq., table value 1.269+0.591i)
  // are kept exactly as given - these are not resonance widths, and
  // have no lineshape-machinery equivalent to substitute.
  //
  // The three G(Q^2) form factors (G_{+000}, G_{+-+0}, G^omega_{+-+0})
  // are the hard-coded Bondar Appendix-B tables from
  // Novosibirsk4pi_GTables.H - genuinely tabulated, not analytic, so
  // there is no lineshape substitute for these either.
  //
  // Bose symmetrization: handled entirely by the explicit permutation
  // sums in Eqs.(4pi-J000-rho/sigma) and (4pi-Jmix-rho/sigma/omega) -
  // these sums are written such that ANY consistent assignment of
  // "which identical pion is q1 vs q3" (mixed-charge mode) or "which
  // is q2 vs q3 vs q4" (neutral mode) gives the same physical current,
  // so FixMode() below does not need any special identical-particle
  // handling beyond identifying the charge/flavour pattern itself.
  //
  // TAUOLA convention note: the source states the published TAUOLA
  // implementation treats the a1pi and omegapi pieces of the mixed-
  // charge mode as INCOHERENT (added at the rate level, not the
  // amplitude level), and that combining them coherently is a model
  // change that "should be labelled as such". This class combines
  // them COHERENTLY (a single summed amplitude, like every other
  // current in this codebase) - a deliberate, flagged departure from
  // the original TAUOLA convention for implementation simplicity.
  class VA_0_Novosibirsk4Pi : public Current_Base {
  private:
    enum class Mode { PiPi0Pi0Pi0, PiPiPiPi0, Unknown };
    Mode m_mode;

    double  m_Mrho, m_Ma1, m_Msigma, m_Momega;
    Total_Width_Base * p_rho_width, * p_a1_width, * p_sigma_width, * p_omega_width;
    double  m_Lambda;
    Complex m_z;
    Complex m_kappaF0;
    double  m_norm;

    Complex DR(const double & s,Total_Width_Base * width,const double & M) const;
    // rho(770) tower PLUS an f0(500)/sigma admixture (Sec. request:
    // "add a set of new models... plus the f0"). Bondar's own t1/t3
    // only ever use a bare D_rho with no admixture at all - this adds
    // one, reusing Sherpa's own registered f0(600) running width, the
    // same "reuse lineshape machinery" choice as everywhere else in
    // this class. NOT from the original paper; m_kappaF0's default is
    // a guess (see SetModelParameters), retune against data.
    Complex DRrho(const double & s) const;
    double  Fa1sq(const double & s) const;
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
