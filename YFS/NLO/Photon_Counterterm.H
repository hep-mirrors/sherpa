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
      Recola only: light-quark masses as Recola has them (0 = massless) and
      its IR scale mu_IR (RECOLA_IR_SCALE).
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

    per external photon (OpenLoops' loopfactor for its pdg 2002 photon).
    On return, log holds a one-line breakdown.
  */
  double OnShellPhotonCounterterm(const Photon_CT_Input &in,
                                  photon_ct_scheme scheme, std::string &log);

  /*!
    Recola versus OpenLoops counterterm, G_mu scheme. They differ only in the
    light-fermion part of dZ_AA = -Pi'_AA(0), so

      c_Recola = c_OL + (alpha_in/4pi) Re{ -Pi_AA^light(M_Z^2)/M_Z^2
                   - Delta alpha - (4/3) sum_f N_c Q_f^2 ln(m_{f,reg}^2/mu^2) } ,

    f over the light fermions (e, mu, tau, u, d, s, c, b), m_{f,reg} = m_f or
    mu_IR for a massless one, mu any common scale (it cancels). The result
    depends on mu_IR (RECOLA_IR_SCALE) through the photon wave function.
  */
  double RecolaPhotonCounterterm(const Photon_CT_Input &in,
                                 photon_ct_scheme scheme, std::string &log);

  /*!
    The imaginary part of the same conversion, needed only at the amplitude
    level (per-helicity real-virtual of CEEX). In the complex-mass scheme it is
    a helicity-independent constant phase from Im dZe(G_mu) and the heavy
    (W, top) part of dZ_AA. Returns 2 Im c_amp/2; 0 for alpha(0). Provider
    independent.
  */
  double PhotonCountertermImaginary(const Photon_CT_Input &in,
                                    photon_ct_scheme scheme);

}

#endif
