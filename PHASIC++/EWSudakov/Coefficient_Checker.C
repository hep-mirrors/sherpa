#include "PHASIC++/EWSudakov/Coefficient_Checker.H"

using namespace PHASIC;
using namespace ATOOLS;

bool Coefficient_Checker::CheckCoeffs(
    const Coeff_Map& coeffs,
    const METOOLS::Spin_Amplitudes& spinampls,
    const Mandelstam_Variables& mandelstam)
{
  auto res = true;
  const auto& refs = ReferenceCoeffs(mandelstam);
  for (const auto& refkv : refs) {
    const auto& type = refkv.first.first;
    if (activecoeffs.find(type) == activecoeffs.end())
      continue;
    const auto& key = refkv.first;
    msg_Debugging() << "Tests for " << key << " reference values:\n";
    for (const auto& helrefpair : refkv.second) {
      const auto& helicities = helrefpair.first;
      const auto idx = spinampls.GetNumber(helicities);
      const auto coeffsit = coeffs.find(key);
      if (coeffsit == coeffs.end())
        THROW(fatal_error, "EW Sudakov coeffs not found");
      if (!CheckCoeff(coeffsit->second[idx], helrefpair.second, helicities))
        res = false;
    }
  }
  return res;
}

bool Coefficient_Checker::CheckCoeff(const Coeff_Value& coeffpair,
                                     Complex ref,
                                     const std::vector<int>& helicities) const
{
  auto res = true;
  for (int i{ 0 }; i < 2; ++i) {
    auto coeff = (i == 0) ? coeffpair.first : coeffpair.second;
    const auto prec = (std::abs(ref) < 10.0) ? 1.e-2 : 1.e-1;
    const auto singlecoeffres
      = (IsBad(coeff.real())
         || std::abs(coeff.real() - ref) < prec);
    if (singlecoeffres) {
      msg_Debugging() << om::green;
    } else {
      msg_Debugging() << om::red;
    }
    if (i == 0) {
      for (const auto& h : helicities)
        msg_Debugging() << h << " ";
      msg_Debugging() << " coeff: " << coeff;
      if (!singlecoeffres)
        res = false;
    }
    if (i == 1) {
      msg_Debugging() << " alt: " << coeff;
    }
    msg_Debugging() << om::reset;
  }
  msg_Debugging() << "\t vs \t  reference value: " << ref << std::endl;
  return res;
}

std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap>
Coefficient_Checker::ReferenceCoeffs(const Mandelstam_Variables& mandelstam)
{
  std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap> coeffs;

  const double u_over_t = mandelstam.u/mandelstam.t;
  const double u_over_s = mandelstam.u/mandelstam.s;
  const double t_over_s = mandelstam.t/mandelstam.s;

  if (procname == "2_2__e-__e+__mu-__mu+") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -2.58;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -4.96;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.96;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.35;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.29;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.37;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.37;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.45;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] = -2.58;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] =  2.58;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] = -1.29;  // -2*R_lq(RL)=-R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] =  1.29;  // +2*R_lq(RL)=+R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] =  1.29;  // same as for RL
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] = -1.29;  // same as for RL
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] =  1.29;  // same as for RL
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] = -9.83;  // (-4*R_lq(LL)-1/(R_lq(LL)*sw^4)/)2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] =  2.88;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] =  2.88;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 0}] = 7.73;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 0}] = 14.9;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 1}] = 14.9;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 1}] = 22.1;

  } else if (procname == "2_2__e-__e+__u__ub") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.86;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -4.25;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.21;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.29;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.50;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.58;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] =  1.72;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] = -1.72;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] =  1.72;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] = -1.72;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] =  0.86;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] = -0.86;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] =  0.86;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] = -0.86;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] =  2.45;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] =  2.45;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] = -10.6;  // -(4*R_lq(LL)+1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] = -10.6;  // -(4*R_lq(LL)+1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 0}] = 5.58;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 0}] = 12.7;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 1}] = 14.0;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 1}] = 21.2;

  } else if (procname == "2_2__e-__e+__t__tb") {
    // same as 2_2__e-__e+__u__ub except for the lYuk coefficients

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.86;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -4.25;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.21;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.29;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.50;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.58;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] =  1.72;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] = -1.72;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] =  1.72;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] = -1.72;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] =  0.86;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] = -0.86;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] =  0.86;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] = -0.86;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] =  2.45;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] =  2.45;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] = -10.6;  // -(4*R_lq(LL)+1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] = -10.6;  // -(4*R_lq(LL)+1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 0}] = 5.58;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 0}] = 12.7;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 1}] = 14.0;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 1}] = 21.2;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 0}] = -10.6;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 0}] = -10.6;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 1}] = -5.30;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 1}] = -5.30;

  } else if (procname == "2_2__e-__e+__d__db") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.43;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -3.82;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.16;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.24;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.67;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.75;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] = -0.86;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] =  0.86;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] = -0.86;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] =  0.86;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] = -0.43;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] =  0.43;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] = -0.43;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] =  0.43;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] = -11.9;  // -(4*R_lq(LL)-1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] = -11.9;  // -(4*R_lq(LL)-1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] =  2.02;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] =  2.02;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 0}] = 4.29;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 0}] = 11.5;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 1}] = 14.0;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 1}] = 21.2;

  } else if (procname == "2_2__e-__e+__b__bb") {
    // same as 2_2__e-__e+__d__db except for the lYuk coefficients

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 0}] = -1.43;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 0}] = -3.82;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 1}] = -4.68;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 1}] = -7.07;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 0}] = 0.16;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 0}] = 0.24;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 1}] = 0.67;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 1}] = 0.75;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 0}] = -0.86;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 0}] =  0.86;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 0}] = -0.86;  // -2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 0}] =  0.86;  // +2*R_lq(RR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 0}] = -0.43;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 0}] =  0.43;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 0}] = -0.43;  // -2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 0}] =  0.43;  // +2*R_lq(RL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 1}] =  0.43;  // -2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 1}] = -0.43;  // +2*R_lq(LR)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 1}] = -11.9;  // -(4*R_lq(LL)-1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 1}] = -11.9;  // -(4*R_lq(LL)-1/(R_lq(LL)*sw^4))/2
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 1}] =  2.02;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 1}] =  2.02;  // -2*R_lq(LL)
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 0}] = 4.29;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 0}] = 11.5;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 1}] = 14.0;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 1}] = 21.2;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 1}] = -5.30;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 1}] = -5.30;

  } else if (procname == "2_2__e-__e+__W+__W-") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 2, 2}] = -4.96;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 2, 2}] = -7.35;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -12.6;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -12.6;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 2, 2}] = 0.37;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 2, 2}] = 0.45;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 1.98;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 1.98;

    // NOTE: t-ch in Sherpa corresponds to u-ch in the Denner/Pozzorini
    // reference (and vice versa), because their process is ordered differently
    // NOTE: if two contributions are given separately, the first is the N-loop
    // and the second the W-loop contribution
    // LT t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] =  4.47;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] =  4.47;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] =  4.47;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] =  4.47;
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = -4.47 - 4.47 * (1.0 - u_over_t);

  } else if (procname == "2_2__e-__e+__P__P") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -8.15;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -8.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 1}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 0}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 0.22;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 0.22;

    // LT t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = 4.47 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] = 4.47 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = 4.47 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] = 4.47 * u_over_s;
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 4.47 * t_over_s;

  } else if (procname == "2_2__e-__e+__Z__P") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -12.2;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -12.2;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 1}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 0}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 0.22;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 0.22;

    // NOTE: 0<->1 and 2<->3 wrt to the Denner/Pozzorini reference, due to a
    // different process ordering
    // LT t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] = 4.47 * (-1.81*t_over_s + u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = 12.56 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] = 4.47 * (-1.81*t_over_s + u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = 12.56 * u_over_s;
    // LT u-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 4.47 * (-1.81*u_over_s + t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 12.56 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 4.47 * (-1.81*u_over_s + t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 12.56 * t_over_s;

  } else if (procname == "2_2__e-__e+__Z__Z") {

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 0, 1}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{0, 0, 1, 0}] = -1.29;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -16.2;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -16.2;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 0, 1}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{0, 0, 1, 0}] = 0.15;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 0.22;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 0.22;

    // LT t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = 12.56 * (u_over_s - 1.81 * t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] = 12.56 * (u_over_s - 1.81 * t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = 12.56 * (u_over_s - 1.81 * t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] = 12.56 * (u_over_s - 1.81 * t_over_s);
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 12.56 * (t_over_s - 1.81 * u_over_s);

  } else {
    THROW(not_implemented, "No test for this proc");
  }
  return coeffs;
}
