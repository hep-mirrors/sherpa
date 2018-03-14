#include "PHASIC++/EWSudakov/Coefficient_Checker.H"

using namespace PHASIC;
using namespace ATOOLS;

bool Coefficient_Checker::CheckCoeffs(
    const Coeff_Map& coeffs,
    const METOOLS::Spin_Amplitudes& spinampls)
{
  auto res = true;
  const auto& refs = ReferenceCoeffs();
  for (const auto& refkv : refs) {
    const auto& type = refkv.first.first;
    if (activecoeffs.find(type) == activecoeffs.end())
      continue;
    const auto& key = refkv.first;
    for (const auto& helrefpair : refkv.second) {
      const auto& helicities = helrefpair.first;
      const auto idx = spinampls.GetNumber(helicities);
      const auto coeffsit = coeffs.find(key);
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
        << key << " coeff: " << coeff
        << "\t vs \t  reference value: " << helrefpair.second
        << om::reset << std::endl;
      if (!singlecoeffres)
        res = false;
    }
  }
  return res;
}

const std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap>&
Coefficient_Checker::ReferenceCoeffs()
{
  static std::map<std::string, std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap>> coeffs;
  if (coeffs.empty()) {
    auto& mapmm = coeffs["2_2__e-__e+__mu-__mu+"];
    mapmm[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -2.58;
    mapmm[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -4.96;
    mapmm[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.96;
    mapmm[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.35;
    mapmm[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.29;
    mapmm[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.37;
    mapmm[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.37;
    mapmm[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.45;
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] =  1.29;  // same as for RL
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] =  1.29;  // same as for RL
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4)
    mapmm[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] =  2.88;  // -4*R_lq(LL)
    mapmm[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] =  2.88;  // -4*R_lq(LL)
    auto& mapuu = coeffs["2_2__e-__e+__u__ub"];
    mapuu[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.86;
    mapuu[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -4.25;
    mapuu[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    mapuu[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    mapuu[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.21;
    mapuu[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.29;
    mapuu[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.50;
    mapuu[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.58;
    auto& mapdd = coeffs["2_2__e-__e+__d__db"];
    mapdd[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.43;
    mapdd[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -3.82;
    mapdd[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    mapdd[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    mapdd[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.16;
    mapdd[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.24;
    mapdd[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.67;
    mapdd[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.75;
    auto& mapWW = coeffs["2_2__e-__e+__W+__W-"];
    mapWW[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 2, 2}] = -4.96;
    mapWW[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 2, 2}] = -7.35;
    mapWW[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -12.6;
    mapWW[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -12.6;
    mapWW[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 2, 2}] = 0.37;
    mapWW[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 2, 2}] = 0.45;
    mapWW[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 1.98;
    mapWW[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 1.98;
    auto& mapPP = coeffs["2_2__e-__e+__P__P"];
    mapPP[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    mapPP[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    mapPP[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -8.15;
    mapPP[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -8.15;
    mapPP[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 1}] = 0.15;
    mapPP[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 0}] = 0.15;
    mapPP[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 0.22;
    mapPP[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 0.22;
    auto& mapZP = coeffs["2_2__e-__e+__Z__P"];
    mapZP[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    mapZP[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    mapZP[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -12.2;
    mapZP[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -12.2;
    mapZP[{EWSudakov_Log_Type::lZ, {}}] = mapPP[{EWSudakov_Log_Type::lZ, {}}];
    auto& mapZZ = coeffs["2_2__e-__e+__Z__Z"];
    mapZZ[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    mapZZ[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    mapZZ[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -16.2;
    mapZZ[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -16.2;
    mapZZ[{EWSudakov_Log_Type::lZ, {}}] = mapPP[{EWSudakov_Log_Type::lZ, {}}];
  }

  // check proc name is inside the few we have
  const auto it = coeffs.find(procname);
  if (it == coeffs.end())
    THROW(not_implemented, "No test for this proc");
  return it->second;
}
