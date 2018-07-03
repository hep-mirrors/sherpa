#ifndef PHASIC__EWSudakov__EWGoupConstants_H
#define PHASIC__EWSudakov__EWGoupConstants_H

#include "PHASIC++/EWSudakov/EWSudakov.H"

namespace PHASIC {

  class EWGroupConstants {

  public:

    EWGroupConstants();

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10, B.25)
    double DiagonalCew(const ATOOLS::Flavour&, int pol) const;

    /// the same, but non-diagonal neutral gauge bosons terms, eq. (B.25)
    double NondiagonalCew() const;

    /// calculate the squared coupling to Z
    double IZ2(const ATOOLS::Flavour&, int pol) const;

    /// calculate the coupling to Z
    double IZ(const ATOOLS::Flavour&, int pol) const;

    /// calculate the coupling to W plus/minus
    Couplings Ipm(const ATOOLS::Flavour& flav, int pol, bool isplus) const;

    /// calculate the one-loop coefficient of the beta function
    double DiagonalBew(const ATOOLS::Flavour& flav, int pol) const;

    /// the same, but non-diagonal neutral gauge bosons terms, eq. (B.41, B.42)
    double NondiagonalBew() const;

    /// sine and cosine (squared) of the Weinberg angle theta_W
    const double m_sw2, m_cw2, m_sw, m_cw;

  };

}

#endif
