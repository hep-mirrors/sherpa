#ifndef MODEL_SM_FourPion_Current_H
#define MODEL_SM_FourPion_Current_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <array>

namespace METOOLS {

  /*!
    The gamma* -> 4 pi current of arXiv:hep-ph/0008262 as updated in
    arXiv:0804.0359 (a1 pi, f0(1370) rho, "rho rho" and omega pi0 terms;
    omega -> 3 pi couplings of arXiv:hep-ph/0512180), written as
      J^mu = sum_i G_i p_i^mu (+ a Q^mu term, dropped: it vanishes against
             every lepton current it meets),
    G_i scalar functions of the pion momenta and of the photon virtuality q2,
    which enters only through the Q^2 towers B_X(q2).

    J_00 (pi+ pi- 2pi0): arguments p[0], p[1] = pi0 pi0, p[2] = pi+, p[3] = pi-.
    J_cc (2pi+ 2pi-), by isospin:
      J_cc(q1,q2,q3,q4) = J_00(q1,q2,q3,q4) + J_00(q4,q2,q3,q1)
                        + J_00(q1,q3,q2,q4) + J_00(q4,q3,q2,q1),
    q2, q3 one charge and q1, q4 the other. The omega term cancels in J_cc.

    Settings: FOUR_PION_FORM_FACTOR, defaults the fit of arXiv:0804.0359:
      a1, rho, rho', rho'', omega, f0: {Mass, Width}; m_pi;
      Coupling_a1, Coupling_f0, Coupling_rhorho, Coupling_omega;
      Tower: {1, 2, 3: {Mass, Width}} the rho', rho'', rho''' of the Q^2 towers;
      Tower_a1, Tower_f0, Tower_omega: their three weights each;
      g_omegarhopi, g_rhopipi, g_rhorho, beta_rho.
  */
  struct FourPion_Parameters {
    double a1m, a1g, rhom, rhog, rho1m, rho1g, rho2m, rho2g, omm, omg, f0m, f0g, pim;
    double c_a1, c_f0, c_rr, c_om;                  // coupl1, coupl2, coupl3, coupl_om
    double tm[3], tg[3];                            // rho', rho'', rho''' of the Q^2 towers
    double b_a1[3], b_f0[3], b_om[3];
    double g_omrhopi, g_rhopipi, g10, beta_t;
    double con;                                     // a1 running-width normalisation
    void Read();
  private:
    static void ReadMW(ATOOLS::Scoped_Settings s,double &m,double &g);
    static void ReadB(ATOOLS::Scoped_Settings s,double b[3]);
  };

  class FourPionCurrent {
  public:
    //! the run's four-pion current (built on first use from the settings)
    static const FourPionCurrent &Instance();

    //! G_i of J_00; omega = false leaves out the omega pi0 term
    std::array<Complex,4> Coefficients00(double q2,const ATOOLS::Vec4D p[4],bool omega) const;
    //! G_i of J_cc, p in the order of the J_cc arguments above
    std::array<Complex,4> CoefficientsCC(double q2,const ATOOLS::Vec4D p[4]) const;
  private:
    FourPion_Parameters m_par;

    FourPionCurrent();
    //! the P-wave rho Breit-Wigner m^2/(m^2 - s - i sqrt(s) Gamma(s)) of the current
    Complex BWrho(double m,double g,double s) const;
    //! [BW(rho) + sum_k b_k BW(rho_k)]/(1 + sum b_k), P-wave widths in Q^2
    Complex Tower(const double b[3],double q2) const;
    //! the omega pi0 term's Q^2 tower times its coupling, c_om B_om(q2)
    Complex OmegaPi(double q2) const;
    //! g_{omega rho pi} g_{rho pi pi} H(s+-, s+0, s-0), H the omega -> rho pi -> 3 pi sum
    Complex OmegaThreePion(double spm,double sp0,double sm0) const;
    //! the a1 -> 3 pi running width shape of arXiv:hep-ph/0008262 (polynomial fit)
    static double GA1(double s,const FourPion_Parameters &p);
  };

}

#endif
