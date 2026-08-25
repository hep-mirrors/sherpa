#ifndef METOOLS_HadronCurrents_RChL_Functions_H
#define METOOLS_HadronCurrents_RChL_Functions_H

#include "ATOOLS/Math/MyComplex.H"

///////////////////////////////////////////////////////////////////////////
//
// Shared analytic building blocks for the Resonance Chiral Lagrangian
// (RChL) hadronic currents of:
//   [1203.3955] O. Shekhovtsova, T. Przedzinski, P. Roig, Z. Was,
//               Phys.Rev. D86 (2012) 113008
//   [1310.1053] I.M. Nugent, T. Przedzinski, P. Roig, O. Shekhovtsova,
//               Z. Was (sigma-meson extension of the 3pi current)
//
// Kept in one place (rather than duplicated in FF_0_PP.C/FF_0_PPP.C)
// because the same H(x,y), A^R/B^R/A^RR/B^RR (axial) and
// C^R/C^RR/D^R/E^R (odd-intrinsic-parity) building blocks of Appendix A,
// Eq.(44) of [1203.3955] appear in both the 3pi and KKpi currents.
//
// KNOWN LIMITATION (flagged rather than silently approximated):
// The two-meson loop function A_PQ(s), Eq.(45)-(49) of [1203.3955],
// needs J_bar'_PQ(0), the derivative at s=0 of J_bar_PQ(s). Neither
// paper spells this constant out in the excerpt we have to work from
// (it says only "coded in file funct_rpt.f", which is not part of the
// text supplied). We set J_tilde_PQ(s) = J_bar_PQ(s) (i.e. drop the
// -s*J_bar'(0) term), which is exact at s=0 and reasonably close near
// threshold, but WILL be numerically off at large s. Anyone doing a
// precision fit to the two-meson vector form factor should not trust
// this piece without cross-checking against funct_rpt.f or re-deriving
// J_bar'(0) analytically.
//
///////////////////////////////////////////////////////////////////////////

namespace METOOLS {
  namespace RChL {

    /////////////////////////////////////////////////////////////////////
    // Eq.(7)-(8) of [1203.3955]: H(x,y) = -lambda0*y + lambda'*x + lambda'',
    // with lambda0 = (lambda' + lambda'')/4.
    /////////////////////////////////////////////////////////////////////
    double H(const double & x,const double & y,
	    const double & lambda_p,const double & lambda_pp);

    // Eq.(8) of [1203.3955]:
    //   lambda'  = F^2 / (2 sqrt(2) FA GV)
    //   lambda'' = -(1 - 2 GV^2 / F^2) lambda'
    // NOTE the second denominator is F, the chiral decay constant - NOT FV.
    // This used to take FV, which with the Table-4 values (F = 0.091337,
    // FV = 0.168652) made lambda'' a factor 2.00 too large and hence both
    // (lambda'+lambda'') and lambda0 = (lambda'+lambda'')/4 a factor 3.4 too
    // SMALL. Those two drive the double-resonance term F1^RR, so the 3pi
    // shape came out badly tilted while the integrated rate stayed close.
    double Lambda_p(const double & F,const double & FA,const double & GV);
    double Lambda_pp(const double & F,const double & GV,
		     const double & lambda_p);

    /////////////////////////////////////////////////////////////////////
    // alpha_2(q2,s1,s2), Eq.(10) of [1203.3955] (pseudoscalar F4 of the
    // 3pi current). Complex because it contains a rho propagator.
    /////////////////////////////////////////////////////////////////////
    Complex alpha2(const double & q2,const double & s1,const double & s2,
		  const double & s3,const double & mpi2,
		  const double & Mrho,const double & GammaRho_s1,
		  const double & GV,const double & FV);

    /////////////////////////////////////////////////////////////////////
    // Appendix A, Eq.(44) of [1203.3955]: axial building blocks for the
    // 3-hadron currents (used by both 3pi and KKpi channels).
    // "GV_FV_ratio" is 2*GV/FV, i.e. the combination (1 - 2GV/FV) that
    // recurs throughout.
    /////////////////////////////////////////////////////////////////////
    double AR(const double & q2,const double & x,const double & y,
	     const double & m1_2,const double & m2_2,const double & m3_2,
	     const double & GV,const double & FV);
    double BR(const double & x,const double & y,
	     const double & m1_2,const double & m2_2,
	     const double & GV,const double & FV);
    double ARR(const double & q2,const double & x,const double & y,
	      const double & m1_2,const double & m2_2,const double & m3_2,
	      const double & lambda_p,const double & lambda_pp);
    double BRR(const double & q2,const double & x,const double & y,
	      const double & z,const double & m1_2,const double & m2_2,
	      const double & m3_2,
	      const double & lambda_p,const double & lambda_pp);

    /////////////////////////////////////////////////////////////////////
    // Appendix A, Eq.(44) of [1203.3955]: odd-intrinsic-parity building
    // blocks for the vector (F5) form factor of KKpi.
    // IMPORTANT: as Eq.(44) itself makes clear, CR/CRR/DR/ER only ever
    // depend on the specific COMBINATIONS of c1..c6/d1..d3/g1..g5 that
    // Table 5 of [1203.3955] actually tabulates (c1-c2+c5, etc.), never
    // on the individual couplings - which is fortunate, since the
    // individual couplings are NOT uniquely fixed by the paper (Table 5
    // only gives combinations). So these take the combinations directly:
    //   c125   = c1-c2+c5
    //   c1256  = c1-c2-c5+2c6
    //   c1238  = c1+c2+8c3-c5
    //   d123   = d1+8d2-d3
    //   g13    = g1-g3
    /////////////////////////////////////////////////////////////////////
    double CR(const double & q2,const double & x,
	     const double & m1_2,const double & m2_2,const double & m3_2,
	     const double & c125,const double & c1256,const double & c1238,
	     const double & c4);
    double CRR(const double & q2,const double & x,const double & m2,
	      const double & d3,const double & d123);
    double DR(const double & q2,const double & x,const double & y,
	     const double & mK2,const double & mpi2,
	     const double & g13,const double & g2,
	     const double & g4,const double & g5);
    double ER(const double & x,const double & y,
	     const double & g13,const double & g2);

    /////////////////////////////////////////////////////////////////////
    // sigma-meson helper functions, Eq.(5)-(7) of [1310.1053].
    /////////////////////////////////////////////////////////////////////
    Complex BWsigma(const double & x,const double & Msigma2,
		    const double & Msigma,const double & Gamma_x);
    // Gamma_sigma(x) = Gamma_sigma * sigma_pi(x)/sigma_pi(Msigma2)
    double  GammaSigma(const double & x,const double & Msigma2,
		       const double & GammaSigma0,const double & mpi2);
    double  Fsigma(const double & q2,const double & x,const double & mpi2,
		   const double & Rsigma2);
    // sigma_P(q2) = sqrt(1-4mP2/q2); returns 0 below threshold instead of
    // going complex, since it only ever appears as a real damping factor
    // here (both in Gamma_sigma and Fsigma's normalisation).
    double  SigmaPhaseSpace(const double & q2,const double & mP2);

    /////////////////////////////////////////////////////////////////////
    // Two-meson loop function A_PQ(s), Eq.(45)-(49) of [1203.3955].
    // See the KNOWN LIMITATION note above the file's header comment:
    // the -s*Jbar'_PQ(0) term is dropped (set to 0). Returns the full
    // complex value; callers wanting Eq.(21)'s "ReA_PQ(s)" should take
    // .real() of the result.
    /////////////////////////////////////////////////////////////////////
    Complex A_PQ(const double & s,const double & mP2,const double & mQ2,
		const double & F,const double & mu2);

    /////////////////////////////////////////////////////////////////////
    // Polynomial fit to the (dominant, pion-channel) q^2-dependent a1
    // width, Eq.(17) of [1310.1053], with the numerical coefficients
    // quoted there (a..p). This is used INSTEAD OF the full double
    // phase-space integral of Eq.(35)-(36) of [1203.3955] (computing
    // that integral live, self-consistently with the very F1/F2/F3
    // being evaluated, is impractical here). The sub-dominant KKpi
    // absorptive-cut contribution to Gamma_a1 (the last two terms of
    // Eq.13 of [1310.1053]) is NOT included - the papers themselves
    // treat it as a pre-tabulated, not-refit correction, and we don't
    // have that table. Below the 3pi threshold (9mpi^2) this returns 0.
    /////////////////////////////////////////////////////////////////////
    double Gamma_a1_PionFit(const double & q2,const double & Mrho,
			   const double & mpi2);
  }
}

#endif
