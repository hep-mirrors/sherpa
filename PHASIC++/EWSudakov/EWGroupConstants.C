#include "PHASIC++/EWSudakov/EWGroupConstants.H"

#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/KF_Table.H"

#include <cassert>

using namespace ATOOLS;
using namespace PHASIC;

EWGroupConstants::EWGroupConstants():
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  // We use the minus sign in \sin\theta_W to translate between the convention
  // of Denner/Pozzorini and the convention of Sherpa. More specifically, it
  // accounts for the sign difference in the definition of the EW covariant
  // derivative in the Denner/Pozzorini language. This can be found in [Böhm,
  // Denner, Joos: Gauge Theories of the Strong and Electroweak Interaction, 3.
  // neubearb. Aufl., Stuttgart 2001], eq. (4.2.4).  For more information, see
  // https://gitlab.com/ebothmann/sherpa/issues/3. By having this minus sign
  // here, we can use all the constants defined by this class in the language
  // of Denner/Pozzorini.
  m_sw{ -sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) },
  m_aew{ MODEL::s_model->ScalarConstant("alpha_QED")},
  m_mw2{sqr(s_kftable[kf_Wplus]->m_mass)},
  m_mw{s_kftable[kf_Wplus]->m_mass},
  m_mz{s_kftable[kf_Z]->m_mass},
  m_mt{s_kftable[kf_t]->m_mass},
  m_mh0{s_kftable[kf_h0]->m_mass},
  m_cvev{MODEL::s_model->ComplexConstant("cvev").real()}
{
  /// Init running values to their default value
  m_ewpar.m_sw2_r  = m_sw2;
  m_ewpar.m_cw2_r  = m_cw2;
  m_ewpar.m_aew_r  = m_aew;
  m_ewpar.m_mw_r   = m_mw;
  m_ewpar.m_mz_r   = m_mz;
  m_ewpar.m_mt_r   = m_mt;
  m_ewpar.m_mh0_r  = m_mh0;
  m_ewpar.m_cvev_r = m_cvev;
}

double EWGroupConstants::dcw2cw2(const double t2) const
{
  // eq. 5.6
  // at the moment all alpha/4pi are taken out
  return m_sw/m_cw*NondiagonalBew()*log(t2/m_mw2);
}

double EWGroupConstants::dsw2sw2(const double t2) const
{
  // eq. 5.6 combined with sw2 = 1-cw2
  // at the moment all alpha/4pi are taken out
  return -m_cw2/m_sw2*dcw2cw2(t2);
}

double EWGroupConstants::dalphaalpha(const double t2) const
{
  return -DiagonalBew(ATOOLS::Flavour(kf_photon),0)*log(t2/m_mw2) + deltaZem();
}

double EWGroupConstants::deltaZem() const
{
  // quarks + leptons
  double res = 0.0;
  for (size_t i=0;i<17; ++i) {
    if (i==7) i=11;
    auto fli = Flavour(i);
    if(!(fli.IsMassive())) continue;
    double temp = (fli.IsQuark()) ? 3 : 1;
    temp*=fli.Charge()*log(m_mw2/sqr(fli.Mass()));
    res+=2./3.*temp;
  }
  return res;
}

double EWGroupConstants::DiagonalCew(const Flavour& flav, int pol) const
{
  assert(!(flav.IsVector() && pol == 2));
  static const auto CewLefthandedLepton = (1 + 2*m_cw2) / (4*m_sw2*m_cw2);
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
    return 2/m_sw2;
  } else if (flav.IsBoson() && flav.Charge() == 0) {
    if (flav.IsPhoton()) {
      return 2.0;
    } else {
      return 2.0 * m_cw2/m_sw2;
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double EWGroupConstants::NondiagonalCew() const noexcept
{
  return -2.0 * m_cw/m_sw;
}

Couplings EWGroupConstants::IZ2(const Flavour& flav, int pol) const
{
  auto couplings{IZ(flav, pol)};
  for (auto& coupling : couplings) {
    coupling.second *= coupling.second;
  }
  return couplings;
}

Couplings EWGroupConstants::IZ(const Flavour& flav, int pol) const
{
  assert(!(flav.IsVector() && pol == 2));
  const auto sign = (flav.IsAnti() ? -1 : 1);
  static const auto IZLefthandedLepton = (m_sw2 - m_cw2)/(2*m_cw*m_sw);
  const long int signed_kf{flav};
  if (std::abs(signed_kf) == kf_phiplus) {
    // add an extra minus sign here wrt the corresponding lepton coupling,
    // because W+ is the particle, whereas W- is the anti-particle, and
    // they correspond via the Goldstone boson equivalence theorem to the
    // positron (anti-particle) and the electron (particle) respectively;
    // i.e. the roles of the particle/anti-particle swap wrt the
    // correspondence
    return {{signed_kf, -sign * IZLefthandedLepton}};
  } else if (signed_kf == kf_chi) {
    return {{kf_h0, {0.0, -1.0 / (2 * m_cw * m_sw)}}};
  } else if (flav.IsLepton()) {
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return {{signed_kf, sign * m_sw / m_cw}};
    } else {
      if (flav.IsUptype())
        return {{signed_kf, sign / (2 * m_sw * m_cw)}};
      else
        return {{signed_kf, sign * IZLefthandedLepton}};
    }
  } else if (flav.IsQuark()) { // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        return {{signed_kf, -sign * 2 / 3.0 * m_sw / m_cw}};
      else
        return {{signed_kf, sign * 1 / 3.0 * m_sw / m_cw}};
    } else {
      if (flav.IsUptype())
        return {{signed_kf, sign * (3 * m_cw2 - m_sw2) / (6 * m_sw * m_cw)}};
      else
        return {{signed_kf, -sign * (3 * m_cw2 + m_sw2) / (6 * m_sw * m_cw)}};
    }
  } else if (std::abs(signed_kf) == kf_Wplus) {
    return {{signed_kf, sign * m_cw / m_sw}};
  } else if (signed_kf == kf_Z) {
    return {};  // the Z self-coupling is zero
  } else if (signed_kf == kf_photon) {
    return {};  // the Z does not couple to the photon
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
  assert(!(flav.IsVector() && pol == 2));
  const long int signed_kf{flav};
  if (flav.IsFermion()) {
    if (pol == 0)
      return {};
    const auto isfermionplus = flav.IsUptype();
    if (flav.IsAnti() && (isplus == isfermionplus))
      return {{flav.IsoWeakPartner(), -1 / (sqrt(2) * m_sw)}};
    else if (!flav.IsAnti() && (isplus != isfermionplus))
      return {{flav.IsoWeakPartner(), 1 / (sqrt(2) * m_sw)}};
    else
      return {};
  } else if (std::abs(signed_kf) == kf_phiplus) {
    if (isplus != flav.IsAnti()) {
      return {};
    }
    return {
	    {kf_chi, {0, -1.0 / (2.0 * m_sw)}},           // I_\chi^\pm
	    {kf_h0, (isplus ? -1.0 : 1.0) / (2.0 * m_sw)} // I_H^\pm
    };
  } else if (std::abs(signed_kf) == kf_Wplus) {
    // cf. (B.22), (B.26) and (B.27)
    if (isplus != flav.IsAnti()) {
      return {};
    }
    return {{kf_photon, isplus ? -1.0 : 1.0},
            {kf_Z, (isplus ? 1.0 : -1.0) * m_cw / m_sw}};
  } else if (signed_kf == kf_chi) {
    return {{(isplus ? 1.0 : -1.0) * kf_phiplus, {0.0, 1.0 / (2.0 * m_sw)}}};
  } else if (signed_kf == kf_Z) {
    // cf. (B.22), (B.26) and (B.27)
    return {{(isplus ? 1.0 : -1.0) * kf_Wplus,
             (isplus ? 1.0 : -1.0) * m_cw / m_sw}};
  } else if (signed_kf == kf_photon) {
    return {{(isplus ? 1.0 : -1.0) * kf_Wplus, (isplus ? -1.0 : 1.0)}};
  } else {
    MyStrStream s;
    s << "Missing implementation for flavour: " << flav
      << " (pol: " << pol << ')';
    THROW(not_implemented, s.str());
  }
}

double EWGroupConstants::DiagonalBew(const ATOOLS::Flavour& flav, int pol) const
{
  const kf_code kf{flav.Kfcode()};
  if (pol != 2) {
    if (kf == kf_Wplus)
      return 19.0/(6.0*m_sw2);
    else if (kf == kf_photon)
      return -11.0/3.0;
    else if (kf == kf_Z)
      return (19.0 - m_sw2*(38.0 + 22.0*m_sw2)) / (6.0*m_sw2*m_cw2);
  }
  MyStrStream s;
  s << "Missing implementation for flavour: " << flav
    << " (pol: " << pol << ')';
  THROW(not_implemented, s.str());
}

double EWGroupConstants::NondiagonalBew() const noexcept
{
  return -(19.0 + 22.0*m_sw2) / (6.0*m_sw*m_cw);
}

MODEL::EWParameters EWGroupConstants::EvolveEWparameters(const double t2) const
{
  m_ewpar.m_cw2_r = m_cw2*(1. + dcw2cw2(t2));
  m_ewpar.m_sw2_r = m_sw2*(1. + dsw2sw2(t2));
  m_ewpar.m_aew_r = m_aew*(1. + dalphaalpha(t2));
  return m_ewpar;
}
