#include "PHASIC++/EWSudakov/Coefficient_Checker.H"

using namespace PHASIC;
using namespace ATOOLS;

bool Coefficient_Checker::CheckCoeffs(
    const std::map<std::string, std::vector<Complex>>& coeffs,
    const METOOLS::Spin_Amplitudes& spinampls)
{
  auto res = true;
  const auto& refs = ReferenceCoeffs();
  for (const auto& typecoeffspair : refs) {
    const auto& type = typecoeffspair.first;
  for (const auto& helrefpair : typecoeffspair.second) {
    const auto& helicities = helrefpair.first;
    const auto idx = spinampls.GetNumber(helicities);
    const auto coeffsit = coeffs.find(type);
    if (coeffsit == coeffs.end())
      THROW(fatal_error, "EW Sudakov coeffs not found");
    const auto coeff = coeffsit->second[idx];
    const auto prec = (std::abs(helrefpair.second) < 10.0) ? 1.e-2 : 1.e-1;
    const auto singlecoeffres =
      (IsBad(coeff.real()) || std::abs(coeff.real() - helrefpair.second) < prec);
    if (singlecoeffres) {
      msg_Debugging() << om::green;
    } else {
      msg_Debugging() << om::red;
    }
    for (const auto& h : helicities)
      msg_Debugging() << h << " ";
    msg_Debugging()
      << type + " coeff: " << coeff
      << "\t vs \t  reference value: " << helrefpair.second
      << om::reset << std::endl;
    if (!singlecoeffres)
      res = false;
  }
}
return res;
}

bool Coefficient_Checker::CheckAngularCoeffs(
    const std::vector<LegIndizes_Coeff_Map>& coeffs,
    const METOOLS::Spin_Amplitudes& spinampls)
{
  auto res = true;
  const auto& refs = AngularReferenceCoeffs();
  for (const auto& legscoeffspair : refs) {
    const auto& legs = legscoeffspair.first;
    for (const auto& helrefpair : legscoeffspair.second) {
      const auto& helicities = helrefpair.first;
      const auto idx = spinampls.GetNumber(helicities);
      const auto legscoeffsmap = coeffs[idx];
      const auto coeffit = legscoeffsmap.find(legs);
      if (coeffit == legscoeffsmap.end())
        THROW(fatal_error, "EW Sudakov coeffs not found");
      const auto coeff = coeffit->second;
      const auto prec = (std::abs(helrefpair.second) < 10.0) ? 1.e-2 : 1.e-1;
      const auto singlecoeffres =
        (IsBad(coeff.real()) || std::abs(coeff.real() - helrefpair.second) < prec);
      if (singlecoeffres) {
        msg_Debugging() << om::green;
      } else {
        msg_Debugging() << om::red;
      }
      for (const auto& h : helicities)
        msg_Debugging() << h << " ";
      msg_Debugging()
        << "(" << legs[0] << ", " << legs[1] << ") coeff: " << coeff
        << "\t vs \t  reference value: " << helrefpair.second
        << om::reset << std::endl;
      if (!singlecoeffres)
        res = false;
    }
  }
  return res;
}

const std::map<std::string, Coefficient_Checker::HelicityCoeffMap>&
Coefficient_Checker::ReferenceCoeffs()
{
  static std::map<std::string, std::map<std::string, HelicityCoeffMap>> coeffs;
  if (coeffs.empty()) {
    auto& mapmm = coeffs["2_2__e-__e+__mu-__mu+"];
    mapmm["L"][{0, 0, 0, 0}] = -2.58;
    mapmm["L"][{1, 1, 0, 0}] = -4.96;
    mapmm["L"][{0, 0, 1, 1}] = -4.96;
    mapmm["L"][{1, 1, 1, 1}] = -7.35;
    mapmm["lZ"][{0, 0, 0, 0}] = 0.29;
    mapmm["lZ"][{1, 1, 0, 0}] = 0.37;
    mapmm["lZ"][{0, 0, 1, 1}] = 0.37;
    mapmm["lZ"][{1, 1, 1, 1}] = 0.45;
    auto& mapuu = coeffs["2_2__e-__e+__u__ub"];
    mapuu["L"][{0, 0, 0, 0}] = -1.86;
    mapuu["L"][{1, 1, 0, 0}] = -4.25;
    mapuu["L"][{0, 0, 1, 1}] = -4.68;
    mapuu["L"][{1, 1, 1, 1}] = -7.07;
    mapuu["lZ"][{0, 0, 0, 0}] = 0.21;
    mapuu["lZ"][{1, 1, 0, 0}] = 0.29;
    mapuu["lZ"][{0, 0, 1, 1}] = 0.50;
    mapuu["lZ"][{1, 1, 1, 1}] = 0.58;
    auto& mapdd = coeffs["2_2__e-__e+__d__db"];
    mapdd["L"][{0, 0, 0, 0}] = -1.43;
    mapdd["L"][{1, 1, 0, 0}] = -3.82;
    mapdd["L"][{0, 0, 1, 1}] = -4.68;
    mapdd["L"][{1, 1, 1, 1}] = -7.07;
    mapdd["lZ"][{0, 0, 0, 0}] = 0.16;
    mapdd["lZ"][{1, 1, 0, 0}] = 0.24;
    mapdd["lZ"][{0, 0, 1, 1}] = 0.67;
    mapdd["lZ"][{1, 1, 1, 1}] = 0.75;
    auto& mapWW = coeffs["2_2__e-__e+__W+__W-"];
    mapWW["L"][{0, 0, 2, 2}] = -4.96;
    mapWW["L"][{1, 1, 2, 2}] = -7.35;
    mapWW["L"][{1, 1, 0, 1}] = -12.6;
    mapWW["L"][{1, 1, 1, 0}] = -12.6;
    mapWW["lZ"][{0, 0, 2, 2}] = 0.37;
    mapWW["lZ"][{1, 1, 2, 2}] = 0.45;
    mapWW["lZ"][{1, 1, 0, 1}] = 1.98;
    mapWW["lZ"][{1, 1, 1, 0}] = 1.98;
    auto& mapPP = coeffs["2_2__e-__e+__P__P"];
    mapPP["L"][{0, 0, 0, 1}] = -1.29;
    mapPP["L"][{0, 0, 1, 0}] = -1.29;
    mapPP["L"][{1, 1, 0, 1}] = -8.15;
    mapPP["L"][{1, 1, 1, 0}] = -8.15;
    mapPP["lZ"][{0, 0, 0, 1}] = 0.15;
    mapPP["lZ"][{0, 0, 1, 0}] = 0.15;
    mapPP["lZ"][{1, 1, 0, 1}] = 0.22;
    mapPP["lZ"][{1, 1, 1, 0}] = 0.22;
    auto& mapZP = coeffs["2_2__e-__e+__Z__P"];
    mapZP["L"][{0, 0, 0, 1}] = -1.29;
    mapZP["L"][{0, 0, 1, 0}] = -1.29;
    mapZP["L"][{1, 1, 0, 1}] = -12.2;
    mapZP["L"][{1, 1, 1, 0}] = -12.2;
    mapZP["lZ"] = mapPP["lZ"];
    auto& mapZZ = coeffs["2_2__e-__e+__Z__Z"];
    mapZZ["L"][{0, 0, 0, 1}] = -1.29;
    mapZZ["L"][{0, 0, 1, 0}] = -1.29;
    mapZZ["L"][{1, 1, 0, 1}] = -16.2;
    mapZZ["L"][{1, 1, 1, 0}] = -16.2;
    mapZZ["lZ"] = mapPP["lZ"];
  }

  // check proc name is inside the few we have
  const auto it = coeffs.find(procname);
  if (it == coeffs.end())
    THROW(not_implemented, "No test for this proc");
  return it->second;
}

const std::map<Two_Leg_Indizes, Coefficient_Checker::HelicityCoeffMap>&
Coefficient_Checker::AngularReferenceCoeffs()
{
  static std::map<std::string, std::map<Two_Leg_Indizes, HelicityCoeffMap>> coeffs;
  if (coeffs.empty()) {
    auto& mapmm = coeffs["2_2__e-__e+__mu-__mu+"];

    mapmm[{2, 0}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    mapmm[{3, 0}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)
    mapmm[{3, 1}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    mapmm[{2, 1}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)

    mapmm[{2, 0}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    mapmm[{3, 0}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    mapmm[{3, 1}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    mapmm[{2, 1}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    mapmm[{2, 0}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    mapmm[{3, 0}][{0, 0, 1, 1}] =  1.29;  // same as for RL
    mapmm[{3, 1}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    mapmm[{2, 1}][{0, 0, 1, 1}] =  1.29;  // same as for RL

    mapmm[{2, 0}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4)
    mapmm[{3, 1}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4)

    mapmm[{3, 0}][{1, 1, 1, 1}] =  2.88;  // -4*R_lq(LL)
    mapmm[{2, 1}][{1, 1, 1, 1}] =  2.88;  // -4*R_lq(LL)
  }

  // check proc name is inside the few we have
  const auto it = coeffs.find(procname);
  if (it == coeffs.end())
    THROW(not_implemented, "No test for this proc");
  return it->second;
}
