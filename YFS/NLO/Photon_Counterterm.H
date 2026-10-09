#ifndef YFS_NLO_Photon_Counterterm_H
#define YFS_NLO_Photon_Counterterm_H

#include <complex>
#include <string>

namespace YFS {

  /*!
    Input of OnShellPhotonCounterterm: the electroweak parameters as the loop
    provider (OpenLoops through Sherpa's interface) sees them. Masses in GeV,
    widths in GeV (complex-mass scheme when cms is true); fermion masses real,
    0 for a massless fermion; u, d, s, c are massless (OpenLoops' libraries).
  */
  //! the loop provider whose external-photon renormalisation the
  //! counterterm converts (see photon_ct_provider::recola below)
  enum class photon_ct_provider { openloops = 0, recola = 1 };

  struct Photon_CT_Input {
    double MZ = 91.1876, WZ = 0., MW = 80.379, WW = 0., MH = 125., WH = 0.,
           MT = 172.5, WT = 0.;
    double me = 0., mmu = 0., mtau = 0., mb = 0.;
    //! alpha(0), and the alpha handed to the provider (Sherpa's model alpha)
    double alpha0 = 1./137.035999, alpha_in = 1./137.035999;
    bool cms = true;
    photon_ct_provider provider = photon_ct_provider::openloops;
    /*!
      Recola only: the light-quark masses as Recola has them (0 = massless,
      then dimensionally regularised), and Recola's IR scale mu_IR
      (RECOLA_IR_SCALE), which is what a massless charged fermion's
      collinear log of the photon wave-function renormalisation is left
      with after the UV cancellation (Delta_IR = 0).
    */
    double mu = 0., md = 0., ms = 0., mc = 0.;
    double mu_ir = 100.;
  };

  //! input scheme of the charge renormalisation the provider applies to the
  //! emitted photon (OpenLoops' ew_renorm_scheme 0, 1, 2)
  enum class photon_ct_scheme { alpha0 = 0, Gmu = 1, alphamZ = 2 };

  /*!
    The relative O(alpha) correction that turns the one-loop |M_{n+1}|^2 of an
    emitted photon renormalised in the input scheme into that of an on-shell
    photon coupling with alpha(0) (tree x alpha(0)/alpha_in not included):

      c = 2 Re[dZe(alpha(0), n-reg) - dZe(scheme)] alpha_in/(4 pi)

    per external photon - exactly OpenLoops' loopfactor for its pdg 2002
    photon (renormalisation_ew.F90, photon_factors, 2711-2764) times the 2
    of 2 Re(M0^* M1). See Photon_Counterterm.C for the pieces.
    On return, log holds a one-line breakdown.
  */
  double OnShellPhotonCounterterm(const Photon_CT_Input &in,
                                  photon_ct_scheme scheme, std::string &log);

  /*!
    Recola 2 (SM_2.2.3 model file) versus OpenLoops, G_mu scheme. Both
    renormalise the emitted photon's vertex with the SAME dZe(G_mu)
    (ctparameter146: dZee = dZeeZM + dr/2, no photon wave-function content,
    Denner's Delta r), and both attach the photon wave-function factor
    dZ_AA/2 to the external photon. They differ ONLY in the light-fermion
    part of dZ_AA = -Pi'_AA(0): OpenLoops replaces it by its "nreg" form
    (photon_factors uses dZAAEWnreg: Pi_light(M_Z^2)/M_Z^2 + Delta alpha
    from the alpha(0)/alpha_in input), Recola's dZA_QED2 (ctparameter75)
    takes it as it is - the exact massive derivative for every fermion
    with a mass, which is 3 ln(m^2/mu^2) in the units below, and 0 for a
    massless fermion, whose scaleless B0'(0;0,0) is (Delta_UV - Delta_IR +
    ln mu_UV^2/mu_IR^2), i.e. the massive form with m -> mu_IR once the
    UV part has cancelled in the total (Delta_IR = 0). Hence

      c_Recola = c_OpenLoops + (alpha_in/4pi) Re[ dZ_AA^nreg - dZ_AA^Recola ]
               = c_OL + (alpha_in/4pi) Re{ -Pi_AA^light(M_Z^2)/M_Z^2
                   - Delta alpha
                   - (4/3) sum_f N_c Q_f^2 ln(m_{f,reg}^2/mu^2) } ,

    f over the fermions OpenLoops treats as light (e, mu, tau, u, d, s, c,
    b), m_{f,reg} = m_f (Recola's mass) or mu_IR for a massless one, mu any
    common scale (it cancels, checked). The heavy set (top, W) is the same
    exact derivative in both and drops out. The remaining mu_IR dependence
    is Recola's: with massless quarks its (n+1)-point amplitude depends on
    the IR scale beyond the soft-photon poles, through the photon wave
    function, so the counterterm carries the same mu_IR (RECOLA_IR_SCALE).
    Z pole, G_mu, e/mu/tau massive, u/d/s/c/b massless, mu_IR = 100 GeV:
    c_OL = 0.0075211, dc = -0.0101, c_Recola = -0.0025.
    Verified by the soft limit of v_{n+1} - v_B (NOTES-ceex-helicity-
    realvirtual-2026-09-29.md), not fitted to it.
  */
  double RecolaPhotonCounterterm(const Photon_CT_Input &in,
                                 photon_ct_scheme scheme, std::string &log);

  /*!
    The imaginary part of the same conversion, needed only at the AMPLITUDE
    level (the per-helicity real-virtual of CEEX). In the complex-mass
    scheme the counterterms are complex, and the emitted photon's vertex
    factor dZe(G_mu) + dZ_AA/2 carries Im[dZe(G_mu)] ~ Im(Delta r)/2, a
    helicity-independent constant phase the Born point does not have; the
    on-shell alpha(0) factor -(sw/cw) Sigma^AZ_T(0)/M_Z^2 has an imaginary
    part of O(alpha Gamma_W/M_W) instead. Per amplitude the conversion is
      c_amp/2 = -(sw/cw) Sigma^AZ_T(0)/M_Z^2 - dZe(G_mu) - dZ_AA^prov/2 ,
    2 Re c_amp/2 = the counterterms above, and this returns 2 Im c_amp/2
    with dZ_AA^prov's imaginary part from its heavy (W, top) pieces (the
    light-fermion parts are real in both providers). Z pole, G_mu:
    +0.0313 (Recola's Im per amplitude 0.0157, measured on the soft limit
    before this term was added). alpha(0): 0. Provider independent.
  */
  double PhotonCountertermImaginary(const Photon_CT_Input &in,
                                    photon_ct_scheme scheme);

}

#endif
