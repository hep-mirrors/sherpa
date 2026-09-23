#ifndef METOOLS_HadronCurrents_FormFactors_Novosibirsk4pi_GTables_H
#define METOOLS_HadronCurrents_FormFactors_Novosibirsk4pi_GTables_H

#include <complex>

namespace METOOLS {
  typedef std::complex<double> Complex;
  // ------------------------------------------------------------------
  // Lookup tables for the Novosibirsk/CVC four-pion current, taken from
  // TAUOLA's own implementation (tauola-factory/Binp/binp.f, TAUOLA
  // 1.1.6c) rather than from Bondar et al., Comput. Phys. Commun. 146
  // (2002) 139.
  //
  // WHY THE CODE AND NOT THE PAPER.  The paper writes the current as
  // J^mu = G(Q^2) x sum(t_i^mu) and prints the three G tables.  That is
  // the STRUCTURE only.  Everything that sets the scale lives in
  // binp.f and appears nowhere in the paper:
  //
  //   * a per-mode constant  (76.565 / 886.838 / 96.867)
  //   * a per-mode threshold factor sqrt(a*Q - b)
  //   * an overall 1/M_rho^4 / Q
  //   * ZFA1TAB, a further interpolated table, on the neutral mode only
  //   * propagators normalised to be dimensionless, and a
  //     Gounaris-Sakurai dispersive real part in the rho that the
  //     paper's Eq.(18) does not have
  //   * t1 carrying the opposite sign to the paper's Eq.(16)
  //
  // A transcription of the paper alone reproduces every equation and
  // still lands about 1e-5 low, with the ratio between the two charge
  // modes wrong as well - the mode dependence above cannot be absorbed
  // into any common constant.
  //
  // Note also that the paper prints its G tables at 96 points; the code
  // carries 98 on the same grid.  These are the code's 98.
  // ------------------------------------------------------------------

  // G(Q) for the three components.  Common grid in Q [GeV], 0.600 to
  // 1.874, step 0.0131313.  binp.f returns 0 below 0.6 and above 1.777
  // (NOT the end of the grid) rather than clamping.
  extern const double kG_arg[98];
  extern const double kG_a1mix[98];   // a1 pi, pi+ pi- pi+ pi0
  extern const double kG_ommix[98];   // omega pi, pi+ pi- pi+ pi0
  extern const double kG_neut[98];    // a1 pi, pi+ pi0 pi0 pi0

  // ZFA1TAB, applied to the neutral mode only.  Argument is Q^2 [GeV^2].
  extern const double kZFA1_arg[100];
  extern const double kZFA1_val[100];


  // Linear interpolation on a uniformly-gridded table.  Returns 0
  // outside [lo,hi] rather than clamping, matching binp.f.
  double Novo4PiLookup(const double & x,const double * arg,const double * val,
                    const int n,const double & lo,const double & hi);

  // ---- the numerical content of the parameterisation -------------
  // These live here rather than in a Total_Width_Base subclass because
  // Total_Width_Base builds a width by summing declared
  // Partial_Width_Base channels and its operator() is not virtual - it
  // cannot express a tabulated width, and forcing one into it would
  // fight the abstraction rather than use it.  They are the numerical
  // content of one parameterisation, so they sit with its tables.

  // a1 form factor.  NOTE scaleA is Lambda^2, not Lambda - binp.f
  // writes 1 + s/Scale_A with Scale_A = 1.2.
  double Novo4PiFa1Sq(const double & s,const double & Ma1,const double & scaleA);
  // Gounaris-Sakurai subtraction for the rho, binp.f dm(s).  Kept here
  // rather than reusing BreitWigner::ValueGS because binp.f's grouping
  // and normalisation differ; see Novo4Pi_Propagator in Propagator.H.
  double Novo4PiHRho(const double & s);
  double Novo4PiDHRho(const double & s);
  double Novo4PiDM(const double & s,const double & Mrho);
}

#endif
