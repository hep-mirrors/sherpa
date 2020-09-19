#include "PHASIC++/EWSudakov/KFactor_Checker.H"
#include "METOOLS/Main/Spin_Structure.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;

bool KFactor_Checker::CheckKFactor(
    double kfac,
    const Mandelstam_Variables& mandelstam,
    const EWGroupConstants& groupconstants)
{
  std::ofstream logfile;
  if (!logfilename.empty())
    logfile.open(logfilename, std::fstream::out | std::fstream::app);
  auto res = true;
  const auto& ref = ReferenceKFactor(mandelstam, groupconstants);
  msg_Debugging() << "Tests for reference values:\n";
  res = CheckKFactor(kfac, ref);
  if (logfile.is_open()) {
    logfile << mandelstam.s << '\t' << mandelstam.t << '\t'
            << mandelstam.u << '\t' << kfac << '\t'
            << ref << '\n';
  }
  return res;
}

bool KFactor_Checker::CheckKFactor(double kfac, double ref) const
{
  auto prec = 0.01;
  const auto res = (IsBad(kfac) || std::abs(kfac - ref) <= prec);
  if (res) {
    msg_Debugging() << om::green;
  } else {
    msg_Debugging() << om::red;
  }
  msg_Debugging() << " kfac: " << kfac;
  msg_Debugging() << om::reset;
  msg_Debugging() << "\t vs \t  reference value: " << ref << std::endl;
  return res;
}

double KFactor_Checker::ReferenceKFactor(const Mandelstam_Variables& mandelstam,
                                         const EWGroupConstants& groupconstants)
{
  double ref {1.0};
  if (procname == "2_2__u__ub__Z__G") {
    // equation adapted from hep-ph/0408308
    Flavour uflav {kf_u};
    const auto ls = std::log(mandelstam.s/groupconstants.m_mw2);
    const auto lt = std::log(std::abs(mandelstam.t)/mandelstam.s);
    const auto lu = std::log(std::abs(mandelstam.u)/mandelstam.s);
    Complex A_0 {0.0}, A_1 {0.0};
    for (int i {0}; i < 2; i++) {
      A_0 += groupconstants.IZ2(uflav, i);
      const auto T3 = (i == 0) ? 0.0 : 1.0/2.0;
      A_1 -= groupconstants.IZ2(uflav, i) *
        groupconstants.DiagonalCew(uflav, i) *
        (sqr(ls) - 3 * ls) +
        sqrt(groupconstants.IZ2(uflav, i)) *
        groupconstants.m_cw / pow(groupconstants.m_sw, 3) * T3 *
        (sqr(ls) + 2 * (lu + lt) * ls);
    }
    ref = 1.0 + (groupconstants.m_aew / (2.0 * M_PI) * (A_1 / A_0)).real();
  } else if (procname == "2_2__d__ub__W-__G") {
    // equation adapted from arXiv:0708.0476
    Flavour uflav {kf_u};
    const auto ls = std::log(mandelstam.s/groupconstants.m_mw2);
    const auto lt = std::log(std::abs(mandelstam.t)/mandelstam.s);
    const auto lu = std::log(std::abs(mandelstam.u)/mandelstam.s);
    Complex A_1 {0.0};
    const auto C_A = 2.0;
    A_1 -= groupconstants.DiagonalCew(uflav, 1) *
      (sqr(ls) - 3 * ls) +
      C_A / (2.0 * groupconstants.m_sw2) *
      (sqr(ls) + 2 * (lu + lt) * ls);
    ref = 1.0 + (groupconstants.m_aew / (2.0 * M_PI) * A_1).real();
  } else if (procname == "2_2__u__ub__P__G") {
    // equation adapted from hep-ph/0508253 (with the replacement Cew - Qq2 ->
    // Cew such that also photon corrections are taken into account)
    Flavour uflav {kf_u};
    const auto ls = std::log(mandelstam.s/groupconstants.m_mw2);
    const auto lt = std::log(std::abs(mandelstam.t)/mandelstam.s);
    const auto lu = std::log(std::abs(mandelstam.u)/mandelstam.s);
    Complex A_0 {2.0 * sqr(uflav.Charge())};
    Complex A_1 {0.0};
    for (int i {0}; i < 2; i++) {
      const auto T3 = (i == 0) ? 0.0 : 1.0/2.0;
      A_1 -= sqr(uflav.Charge()) *
        //(groupconstants.DiagonalCew(uflav, i) - sqr(uflav.Charge())) *
        groupconstants.DiagonalCew(uflav, i) *
        (sqr(ls) - 3 * ls) +
        uflav.Charge() * T3 / groupconstants.m_sw2 *
        (sqr(ls) + 2 * (lu + lt) * ls);
    }
    ref = 1.0 + (groupconstants.m_aew / (2.0 * M_PI) * (A_1 / A_0)).real();
  } else {
    THROW(not_implemented, "No test for this proc");
  }
  return ref;
}
