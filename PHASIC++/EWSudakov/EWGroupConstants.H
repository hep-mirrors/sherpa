#ifndef PHASIC__EWSudakov__EWGoupConstants_H
#define PHASIC__EWSudakov__EWGoupConstants_H

#include "PHASIC++/EWSudakov/EWSudakov.H"

namespace PHASIC {

  class EWGroupConstants {

  public:

    EWGroupConstants();

    /// Returns a diagonal value of the effective electroweak Casimir operator.
    ///
    /// Cf. eq. (B.10).
    ///
    /// \param flav The particle specification $i_k$ in $C^{ew}_{i_k i_k}(k)$.
    /// Note that longitudinal gauge bosons are not allowed. Use the Goldstone
    /// equivalence theorem eq. (3.4) instead and pass the corresponding
    /// Goldstone boson.
    /// \param pol The chirality/polarisation of the particle:
    /// 0: + (right-handed transverse polarisation)
    /// 1: - (left-handed transverse polarisation)
    /// 2: 0 (longitudinal polarisation)
    /// Note that the polarisation will be ignored in the case of Goldstone
    /// bosons. Hence the user does not need to use the correct polarisation
    /// when the Goldstone equivalence theorem is used.
    double DiagonalCew(const ATOOLS::Flavour& flav, int pol) const;

    /// Returns the value of the non-diagonal (neutral gauge-boson mixing)
    /// effective electroweak Casimir operator, cf. eq. (B.25).
    double NondiagonalCew() const noexcept;

    /// Returns IZ() with the coupling values squared.
    Couplings IZ2(const ATOOLS::Flavour&, int pol) const;

    /// Returns the couplings with the Z boson, cf. eq. (B.7).
    ///
    /// The parameter descriptions of DiagonalCew() apply here, too.
    Couplings IZ(const ATOOLS::Flavour&, int pol) const;

    /// Returns the couplings with the W boson.
    ///
    /// The parameter descriptions of DiagonalCew() apply here, too.
    ///
    /// \param isplus Whether the W boson (outgoing) carries away a positive
    /// change.
    Couplings Ipm(const ATOOLS::Flavour& flav, int pol, bool isplus) const;

    /// Returns a diagonal beta function coefficient, cf. eqs. (B.38-42)
    ///
    /// \param flav The particle specification $i_k$ in $b^{ew}_{i_k i_k}(k)$.
    /// For now, only transverse vector bosons are supported.
    /// \param pol Cf. the description of \p pol of DiagonalCew(), and note
    /// that so far only transverse polarisations are allowed.
    double DiagonalBew(const ATOOLS::Flavour& flav, int pol) const;

    /// Returns a non-diagonal beta function coefficient, cf. eqs. (B.42)
    double NondiagonalBew() const noexcept;

    /// sine and cosine (squared) of the Weinberg angle theta_W
    const double m_sw2, m_cw2, m_sw, m_cw;

  };

}

#endif
