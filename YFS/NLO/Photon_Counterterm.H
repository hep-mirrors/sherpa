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
  struct Photon_CT_Input {
    double MZ = 91.1876, WZ = 0., MW = 80.379, WW = 0., MH = 125., WH = 0.,
           MT = 172.5, WT = 0.;
    double me = 0., mmu = 0., mtau = 0., mb = 0.;
    //! alpha(0), and the alpha handed to the provider (Sherpa's model alpha)
    double alpha0 = 1./137.035999, alpha_in = 1./137.035999;
    bool cms = true;
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

}

#endif
