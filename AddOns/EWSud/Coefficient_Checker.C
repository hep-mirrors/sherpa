#include "AddOns/EWSud/Coefficient_Checker.H"
#include "METOOLS/Main/Spin_Structure.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace EWSud;

bool Coefficient_Checker::CheckCoeffs(
    const Coeff_Map& coeffs,
    const Spin_Amplitudes& spinampls,
    const Mandelstam_Variables& mandelstam,
    const EWGroupConstants& groupconstants)
{
  std::ofstream logfile;
  if (!logfilename.empty())
    logfile.open(logfilename, std::fstream::out | std::fstream::app);
  auto res = true;
  const auto& refs = ReferenceCoeffs(mandelstam, groupconstants);
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
        THROW(fatal_error, "EWSud coeffs not found.");
      if (!CheckCoeff(coeffsit->second[idx], helrefpair.second,
		      helicities, coeffsit->first.first))
        res = false;
      if (logfile.is_open()) {
        logfile << key << '\t';
        for (const auto& h : helicities) {
          logfile << h;
        }
        logfile << '\t' << mandelstam.s << '\t' << mandelstam.t << '\t'
                << mandelstam.u << '\t' << coeffsit->second[idx] << '\t'
                << helrefpair.second << '\n';
      }
    }
  }
  return res;
}

bool Coefficient_Checker::CheckCoeff(const Coeff_Value& coeff, Complex ref,
                                     const std::vector<int>& helicities,
				     const EWSudakov_Log_Type ewlt) const
{
  auto res = true;
  auto prec = std::max(std::abs(ref.real())*0.1, 0.05);
  if (ewlt == EWSudakov_Log_Type::lPR) prec = std::abs(ref.real())*0.3;
  const auto singlecoeffres =
      (IsBad(coeff.real()) || std::abs(coeff.real() - ref) <= prec);
  if (singlecoeffres) {
    msg_Debugging() << om::green;
  } else {
    msg_Debugging() << om::red;
  }
  for (const auto& h : helicities) {
    msg_Debugging() << h << " ";
  }
  msg_Debugging() << " coeff: " << coeff;
  if (!singlecoeffres) {
    res = false;
  }
  msg_Debugging() << om::reset;
  msg_Debugging() << "\t vs \t  reference value: " << ref << std::endl;
  return res;
}

std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap>
Coefficient_Checker::ReferenceCoeffs(const Mandelstam_Variables& mandelstam,
                                     const EWGroupConstants& groupconstants)
{
  std::map<Coeff_Map_Key, Coefficient_Checker::HelicityCoeffMap> coeffs;

  // auxiliary quantities needed to calculate some of the coeffs below
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
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 1}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 1}] = -9.03;

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
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 1}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 1}] = -9.03;

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
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 1}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 1}] = -9.03;

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
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 1}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 1}] = -16.6;

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
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 0}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 1}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 1}] = -16.6;

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
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 2, 2}] =  1.29;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 2, 2}] =  1.29;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 2, 2}] =  2.88;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 2, 2}] =  2.88;
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = -4.47 - 4.47 * (1.0 - u_over_t);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 2, 2}] = -1.29;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 2, 2}] = -1.29;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 2, 2}] = -9.83;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 2, 2}] = -9.83;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 2, 2}] = 18.6;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 2, 2}] = 25.7;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 25.2;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 25.2;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 2, 2}] = -31.8;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 2, 2}] = -31.8;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 2, 2}] = 8.80;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 2, 2}] = -9.03;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = -14.20;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = -14.20;


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
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 0}] = 0.0;
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 4.47 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 0}] = 0.20;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 1}] = 0.20;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 7.36;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 7.36;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 0}] = 3.67;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 1}] = 3.67;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = 3.67;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = 3.67;

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
    // different process ordering; i.e. t-ch and u-ch assignment does not
    // change, but {3, 1} here refers to {2, 0} in Denner/Pozzorini and so on
    // LT t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] = 4.47 * (-1.81*t_over_s + u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = 12.56 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] = 4.47 * (-1.81*t_over_s + u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = 12.56 * u_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 0}] = 0.0;
    // LT u-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 4.47 * (-1.81*u_over_s + t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 12.56 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 4.47 * (-1.81*u_over_s + t_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 12.56 * t_over_s;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 0}] = -11.3;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 1}] = -11.3;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 28.1;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 28.1;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 0}] = 15.10;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 1}] = 15.10;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = -17.10;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = -17.10;
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
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{0, 0, 1, 0}] = 0.0;
    // LT u-ch
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = 12.56 * (t_over_s - 1.81 * u_over_s);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{0, 0, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 1, 0}] = -22.8;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{0, 0, 0, 1}] = -22.8;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 48.9;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 48.9;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{0, 0, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 1, 0}] = 26.60;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{0, 0, 0, 1}] = 26.60;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = -37.9;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = -37.9;

  } else if (procname == "2_2__u__db__Z__W+") {

    // auxiliary quantities needed for the u db -> Z W+ coefficients
    const double F_plus  = 1.0/t_over_s + 1.0/u_over_s;
    const double F_minus = 1.0/t_over_s - 1.0/u_over_s;
    const double sw2 = groupconstants.m_sw2;
    const double cw2 = groupconstants.m_cw2;
    // NOTE: we do not use the sign-flipped sinW here, because the coeffs are
    // evaluated within Denner/Pozzorini's conventions, so we need to be
    // consistent
    const double sw = groupconstants.m_sw;
    const double cw = groupconstants.m_cw;
    const double YqL = 1.0/3.0;
    const double HZ = (cw2 * F_minus - sw2 * YqL * F_plus) / (2 * sw * cw);
    const double HA = - (F_minus + YqL * F_plus)/2.0;
    const double GZ_plus = cw2 * F_plus / (cw2 * F_minus - sw2 * YqL * F_plus);
    const double GZ_minus = cw2 * F_minus / (cw2 * F_minus - sw2 * YqL * F_plus);

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 2, 2}] = -7.07;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -7.86 - 4.47 * GZ_minus;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -7.86 - 4.47 * GZ_minus;

    // t-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = -4.47 - 4.91 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] = -4.47 - 4.04 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = -4.47 - 4.91 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] = -4.47 - 4.04 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 2, 2}] = -4.47 + 0.43;  // NOTE: A,Z contrib is -2.02
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 2, 2}] = -4.47 + 0.43;  // NOTE: A,Z contrib is -2.02

    // u-ch;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = -4.47 + 4.90 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = -4.47 + 4.05 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = -4.47 + 4.90 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = -4.47 + 4.05 * GZ_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 2, 2}] = -4.47 - 0.43;  // NOTE: A,Z contrib is -2.46
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 2, 2}] = -4.47 - 0.43;  // NOTE: A,Z contrib is -2.46

    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 2, 2}] = 0.92;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 1.32;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 1.32;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 2, 2}] = 24.88;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 21.77 - 9.57 * HA/HZ;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 21.77 - 9.57 * HA/HZ;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 2, 2}] = -31.83;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 2, 2}] = -14.16;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = -11.60 + 9.57 * HA/HZ;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = -11.60 + 9.57 * HA/HZ;

  } else if (procname == "2_2__u__db__W+__P") {

    // NOTE: compared to Pozzorini's thesis, the initial-state quarks are
    // flipped, i.e. t- und u-channel are swapped

    // auxiliary quantities needed for the u db -> Z W+ coefficients
    const double F_plus  = 1.0/u_over_s + 1.0/t_over_s;
    const double F_minus = 1.0/u_over_s - 1.0/t_over_s;
    const double sw2 = groupconstants.m_sw2;
    const double cw2 = groupconstants.m_cw2;
    // NOTE: we do not use the sign-flipped sinW here, because the coeffs are
    // evaluated within Denner/Pozzorini's conventions, so we need to be
    // consistent
    const double sw = groupconstants.m_sw;
    const double cw = groupconstants.m_cw;
    const double YqL = 1.0/3.0;
    const double GA_plus =  F_plus  / (F_minus + YqL * F_plus);
    const double GA_minus = F_minus / (F_minus + YqL * F_plus);

    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 0, 1}] = -7.86 - 4.47 * GA_minus;
    coeffs[{EWSudakov_Log_Type::Ls, {}}][{1, 1, 1, 0}] = -7.86 - 4.47 * GA_minus;

    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 0, 1}] = -3.81 - 0.67 - 5.96 * GA_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 1}}][{1, 1, 1, 0}] = -3.81 - 0.67 - 5.96 * GA_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 0, 1}] = -4.47 * (GA_plus + GA_minus);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 0}}][{1, 1, 1, 0}] = -4.47 * (GA_plus + GA_minus);

    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 0, 1}] = -3.14 - 1.33 + 2.98 * GA_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {2, 0}}][{1, 1, 1, 0}] = -3.14 - 1.33 + 2.98 * GA_plus;
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 0, 1}] =  4.47 * (GA_plus - GA_minus);
    coeffs[{EWSudakov_Log_Type::lSSC, {3, 1}}][{1, 1, 1, 0}] =  4.47 * (GA_plus - GA_minus);

    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 0, 1}] = 1.32;
    coeffs[{EWSudakov_Log_Type::lZ, {}}][{1, 1, 1, 0}] = 1.32;

    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 0, 1}] = 15.42;
    coeffs[{EWSudakov_Log_Type::lC, {}}][{1, 1, 1, 0}] = 15.42;

    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 0, 1}] = 0.0;
    coeffs[{EWSudakov_Log_Type::lYuk, {}}][{1, 1, 1, 0}] = 0.0;

    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 0, 1}] = -5.25;
    coeffs[{EWSudakov_Log_Type::lPR, {}}][{1, 1, 1, 0}] = -5.25;

  } else {
    THROW(not_implemented, "No EWSudakov coeff test for this proc: " + procname);
  }
  return coeffs;
}
