#include "PHASIC++/EWSudakov/EWGroupConstants.H"

#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"

using namespace ATOOLS;
using namespace PHASIC;

EWGroupConstants::EWGroupConstants():
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) }
{}

double EWGroupConstants::DiagonalCew(const Flavour& flav, int pol) const
{
  // pol is either chirality or polarisation:
  // 0: + (right-handed or transverse polarisation)
  // 1: - (left-handed or transverse polarisation)
  // 2: 0 (longitudinal polarisation)
  // NOTE: for longitudinal bosons, use the Goldstone equivalence theorem
  static auto CewLefthandedLepton = (1 + 2*m_cw2) / (4*m_sw2*m_cw2);
  if (flav.IsLepton()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return 1/m_cw2;
    } else {
      return CewLefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 1) {
      return (m_sw2 + 27*m_cw2) / (36*m_sw2*m_cw2);
    } else {
      if (flav.IsUptype())
        return 4 / (9*m_cw2);
      else
        return 1 / (9*m_cw2);
    }
  } else if (flav.IsScalar()) {  // cf. eq. (B.18) and (B.16)
    return CewLefthandedLepton;
  } else if (flav.Kfcode() == kf_Wplus) {
    if (pol == 2)
      return CewLefthandedLepton;
    else
      return 2/m_sw2;
  } else if (flav.IsBoson() && flav.Charge() == 0) {
    if (pol == 2) {
      assert(!flav.IsPhoton());
      return CewLefthandedLepton;
    } else if (flav.IsPhoton()) {
      return 2.0;
    } else {
      return 2.0 * m_cw2/m_sw2;
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double EWGroupConstants::NondiagonalCew() const
{
  return -2.0 * m_cw/m_sw;
}

double EWGroupConstants::IZ2(const Flavour& flav, int pol) const
{
  // TODO: just return IZ()^2 here, as soon as it implements the photon/Z
  // constants
  static auto IZ2LefthandedLepton
    = std::pow(m_cw2 - m_sw2, 2) / (4*m_sw2*m_cw2);
  static auto IZ2Neutrino = 1 / (4*m_sw2*m_cw2);
  if (flav.IsLepton()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return m_sw2/m_cw2;
    } else {
      if (flav.IsUptype())
        return IZ2Neutrino;
      else
        return IZ2LefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        return 4*m_sw2 / (9*m_cw2);
      else
        return 1*m_sw2 / (9*m_cw2);
    } else {
      if (flav.IsUptype())
        return std::pow(3*m_cw2 - m_sw2, 2) / (36*m_sw2*m_cw2);
      else
        return std::pow(3*m_cw2 + m_sw2, 2) / (36*m_sw2*m_cw2);
    }
  } else if (flav.Kfcode() == kf_Wplus) {
    if (pol == 2)
      return IZ2LefthandedLepton;
    else
      return m_cw2/m_sw2;
  } else if (flav.IsBoson() && flav.Charge() == 0) {
    if (pol == 2) {
      assert(!flav.IsPhoton());
      return IZ2Neutrino;
    } else {
      return 0.0;
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double EWGroupConstants::IZ(const Flavour& flav, int pol) const
{
  const auto sign = (flav.IsAnti() ? -1 : 1);
  static auto IZLefthandedLepton = (m_sw2 - m_cw2)/(2*m_cw*m_sw);
  if (flav.IsScalar())
    THROW(not_implemented,
          "non-diagonal Z coupling terms for scalars not implemented");
  if (flav.IsLepton()) {
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return sign * m_sw/m_cw;
    } else {
      if (flav.IsUptype())
        return sign / (2*m_sw*m_cw);
      else
        return sign * IZLefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        return -sign * 2/3.0 * m_sw/m_cw;
      else
        return sign * 1/3.0 * m_sw/m_cw;
    } else {
      if (flav.IsUptype())
        return sign * (3*m_cw2 - m_sw2) / (6*m_sw*m_cw);
      else
        return -sign * (3*m_cw2 + m_sw2) / (6*m_sw*m_cw);
    }
  } else if (flav.Kfcode() == kf_Wplus) {
    if (pol == 2) {
      // add an extra minus sign here wrt the corresponding lepton coupling,
      // because W+ is the particle, whereas W- is the anti-particle, and
      // they correspond via the Goldstone boson equivalence theorem to the
      // positron (anti-particle) and the electron (particle) respectively;
      // i.e. the roles of the particle/anti-particle swap wrt the
      // correspondence
      return -sign * IZLefthandedLepton;
    } else {
      return sign * m_cw/m_sw;
    }
  } else if (flav.Kfcode() == kf_Z) {
    if (pol == 2) {

      THROW(not_implemented,
            "non-diagonal Z coupling terms for longitudinal Z not implemented");

      // TODO: enable the following snippet, which requires to modify the
      // signature of this function to return `Couplings`

      // we return the coupling to the scalar, this is corrected by multiplying
      // here with an extra factor of i, cf. (4.26)
      //return {kf_h0, 1.0 / (2.0*m_sw*m_cw)};

    } else {
      return 0.0;  // the Z self-coupling is zero
    }
  } else if (flav.Kfcode() == kf_photon) {
    return 0.0;  // the Z does not couple to the photon
  } else {
    MyStrStream s;
    s << "Missing implementation for flavour: " << flav;
    THROW(not_implemented, s.str());
  }
}

Couplings EWGroupConstants::Ipm(const Flavour& flav,
                                int pol,
                                bool isplus) const
{
  if (flav.IsFermion()) {
    if (pol == 0)
      return {};
    const auto isfermionplus = flav.IsUptype();
    if (flav.IsAnti() && (isplus == isfermionplus))
      return { {flav.IsoWeakPartner().Kfcode(), -1 / (sqrt(2)*m_sw)} };
    else if (!flav.IsAnti() && (isplus != isfermionplus))
      return { {flav.IsoWeakPartner().Kfcode(),  1 / (sqrt(2)*m_sw)} };
    else
      return {};
  } else if (flav.Kfcode() == kf_Wplus) {
    // cf. (B.22), (B.26) and (B.27)
    if (isplus != flav.IsAnti())
      return {};
    if (pol == 2) {
      return {
        // we return the coupling to the pseudoscalar, but tell the recipient
        // to use the ME with the Z instead of the W which makes use of the
        // Goldstone equivalence theorem; this is corrected by multiplying here
        // with an extra factor of (-i), cf. (4.26)
        {kf_Z, -1.0 / (2.0*m_sw)},  // -i * I_\chi^\pm
        {kf_h0, (isplus ? -1.0 : 1.0) / (2.0*m_sw)}  // I_H^\pm
      };
    } else {
      return {
        {kf_photon, isplus ? -1.0 : 1.0},
        {kf_Z, (isplus ? 1.0 : -1.0) * m_cw/m_sw}
      };
    }
  } else if (flav.Kfcode() == kf_Z) {
    // cf. (B.22), (B.26) and (B.27)
    if (pol == 2) {
      // we assume the incoming flavour is the \chi instead of the Z in
      // accordance with the Goldstone equivalence theorem; we correct this
      // by multiplying an extra factor of i, cf. (4.26)
      return { {kf_Wplus, -1.0 / (2.0*m_sw)} };
    } else {
      return { {kf_Wplus, (isplus ? -1.0 : 1.0) * m_cw/m_sw} };
    }
  } else if (flav.Kfcode() == kf_photon) {
    return { {kf_Wplus, (isplus ? 1.0 : -1.0)} };
  } else {
    MyStrStream s;
    s << "Missing implementation for flavour: " << flav
      << " (pol: " << pol << ')';
    THROW(not_implemented, s.str());
  }
}

double EWGroupConstants::DiagonalBew(const ATOOLS::Flavour& flav, int pol) const
{
  if (pol != 2) {
    if (flav.Kfcode() == kf_Wplus)
      return 19.0/(6.0*m_sw2);
    else if (flav.Kfcode() == kf_photon)
      return -11.0/3.0;
    else if (flav.Kfcode() == kf_Z)
      return (19.0 - m_sw2*(38.0 + 22.0*m_sw2)) / (6.0*m_sw2*m_cw2);
  }
  MyStrStream s;
  s << "Missing implementation for flavour: " << flav
    << " (pol: " << pol << ')';
  THROW(not_implemented, s.str());
}

double EWGroupConstants::NondiagonalBew() const
{
  return -(19.0 + 22.0*m_sw2) / (6.0*m_sw*m_cw);
}
