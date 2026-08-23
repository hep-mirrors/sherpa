#ifndef METOOLS_HadronCurrents_FormFactors_Novosibirsk4pi_GTables_H
#define METOOLS_HadronCurrents_FormFactors_Novosibirsk4pi_GTables_H

namespace METOOLS {
  // ------------------------------------------------------------------
  // Novosibirsk/CVC 4pi current lookup tables (Bondar et al., Appendix
  // B, Table 4 - reproduced in tau_two_meson_currents_KS_RChiT.tex,
  // Appendix "Novosibirsk 4pi lookup tables"). These are NOT
  // Breit-Wigner ansaetze - the source note is explicit that a
  // faithful implementation must use these tabulated values verbatim
  // rather than an invented resonance sum (a cross-check on this
  // conversation found the SHAPE of these tables - an abrupt turn-on
  // right at the a1pi/omegapi threshold, peaking immediately, then a
  // smooth monotonic decline with no second local maximum anywhere
  // near where rho'(1450)/rho''(1700) would sit - inconsistent with a
  // simple rho'+rho'' interpretation, so there was never a good
  // shortcut to replace these with two resonances anyway).
  //
  // Three separate tables, each 96 points on the SAME Q grid from
  // 0.600 to 1.874 GeV (Q in GeV, the four-pion invariant mass - the
  // original table mislabels the column Q^2, but the numerical
  // entries are plainly the linear energy grid, not its square):
  //   kG_pi000      : G_{pi+pi0pi0pi0}      (a1pi component, Eq.4pi-J000-rho/sigma)
  //   kG_pimix      : G_{pi+pi-pi+pi0}      (a1pi component, Eq.4pi-Jmix-rho/sigma)
  //   kGomega_pimix : G^omega_{pi+pi-pi+pi0} (omegapi component, Eq.4pi-Jmix-omega)
  //
  // NOTE ON THE HIGH-Q END: all three tables show an abrupt small
  // upward jump at Q=1.769 GeV (most pronounced in kG_pimix: 0.493 ->
  // 0.531 -> 0.857) immediately followed by a hard drop to exactly
  // 0.000 for every remaining grid point (Q=1.782 through 1.874).
  // This pattern is IDENTICAL in shape across all three independently
  // tabulated functions, which makes a transcription error on this
  // end unlikely - it is almost certainly a genuine artifact of the
  // original Bondar et al. table (e.g. the fit becoming unreliable
  // right at the edge of their fitted range before being forced to
  // zero at the kinematic endpoint). Transcribed verbatim regardless,
  // per the source note's own explicit instruction not to alter the
  // tabulated values.
  // ------------------------------------------------------------------

  extern const double kG_pi000[96][2];
  extern const double kGomega_pimix[96][2];
  extern const double kG_pimix[96][2];

  // Linear interpolation lookup on one of the three tables above.
  // Boundary policy (the source note explicitly leaves this to the
  // implementation - "the generator should follow its chosen
  // boundary/interpolation policy explicitly" - rather than silently
  // picking something): CLAMP to the nearest tabulated endpoint value
  // for Q outside [0.600,1.874] GeV, rather than extrapolating. This
  // is a plain, explicit choice, not necessarily identical to
  // whatever the original TAUOLA implementation does at the
  // boundary - flagged here so it is easy to revisit if a different
  // policy is wanted.
  double NovosibirskGLookup(const double & Q, const double table[96][2]);
}

#endif
