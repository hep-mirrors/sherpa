#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"
#include <unordered_map>
#include "SHERPA/SoftPhysics/LDME.H"
#include "ATOOLS/Org/Scoped_Settings.H"

std::unordered_map<kf_code, double> ldme_map;

using namespace ATOOLS;


void msgLDME(kf_code kfc, double LDME)
{
  msg_Info()
    <<"LDME for kf_code "
    << kfc
    << " has been set to "
    << LDME<<"\n";
}

void LoadLDME()
{
  static bool loaded = false;
  if (loaded) return;

  //OCTETS values from tune
  //CCBAR
  auto onia = Settings::GetMainSettings()["QUARKONIA"];
  onia["1S0_441"].SetDefault(0.1E-03);
  onia["1S0_443"].SetDefault(15.081E-03);
  onia["1S0_100443"].SetDefault(21.395E-03);
  onia["3S1_441"].SetDefault(0.1E-03);
  onia["3S1_443"].SetDefault(3.430E-03);
  onia["3S1_100443"].SetDefault(2.628E-03);
  onia["3S1_10441"].SetDefault(2.704E-03);
  double ldme_3S1_c_8_chi_c0_1P = onia["3S1_10441"].Get<double>();
  onia["3S1_20443"].SetDefault(3.*ldme_3S1_c_8_chi_c0_1P);
  onia["3S1_445"].SetDefault(5.*ldme_3S1_c_8_chi_c0_1P);
  onia["3P0_443"].SetDefault(0.0);
  double ldme_3P0_443 = onia["3P0_443"].Get<double>();
  onia["3P0_100443"].SetDefault(0.0);
  double ldme_3P0_100443 = onia["3P0_100443"].Get<double>();
  onia["3P1_443"].SetDefault(3.*ldme_3P0_443);
  onia["3P1_100443"].SetDefault(3.*ldme_3P0_100443);
  onia["3P2_443"].SetDefault(5.*ldme_3P0_443);
  onia["3P2_100443"].SetDefault(5.*ldme_3P0_100443);

  //BBBAR


  //---------------------------------------------------
  //SINGLETS   values from hep-ph/9503356
  //CCBAR
  onia["441"].SetDefault(3.*3.*0.810/(2.*M_PI));
  onia["443"].SetDefault(3.*3.*0.810/(2.*M_PI));
  onia["100443"].SetDefault(3.*3.*0.529/(2.*M_PI));
  onia["10441"].SetDefault(3.*3.*0.075/(2.*M_PI));
  onia["20443"].SetDefault(3.*3.*3.*0.075/(2.*M_PI));
  onia["445"].SetDefault(5.*3.*3.*0.075/(2.*M_PI));
  //BBBAR
  onia["551"].SetDefault(3.*3.*6.477/(2.*M_PI));
  onia["553"].SetDefault(3.*3.*6.477/(2.*M_PI));
  onia["100553"].SetDefault(3.*3.*3.234/(2.*M_PI));
  onia["200553"].SetDefault(3.*3.*3.*2.474/(2.*M_PI));
  onia["10551"].SetDefault(3.*3.*1.417/(2.*M_PI));
  onia["20553"].SetDefault(3.*3.*3.*1.417/(2.*M_PI));
  onia["555"].SetDefault(5.*3.*3.*1.417/(2.*M_PI));
  onia["110551"].SetDefault(3.*3.*1.654/(2.*M_PI));
  onia["120553"].SetDefault(3.*3.*3.*1.654/(2.*M_PI));
  onia["100555"].SetDefault(5.*3.*3.*1.654/(2.*M_PI));


  auto load = [&](auto& group, kf_code kfc, const char* key)
  {
    auto entry = group[key];
    ldme_map[kfc] = entry.template Get<double>();
    if (entry.IsSetExplicitly())
      msgLDME(kfc, ldme_map[kfc]);
  };

  // --- Singlets ---
  load(onia, kf_eta_c_1S,  "441");
  load(onia, kf_J_psi_1S,  "443");
  load(onia, kf_psi_2S,    "100443");
  load(onia, kf_chi_c0_1P, "10441");
  load(onia, kf_chi_c1_1P, "20443");
  load(onia, kf_chi_c2_1P, "445");

  load(onia, kf_eta_b,     "551");
  load(onia, kf_Upsilon_1S,"553");
  load(onia, kf_Upsilon_2S,"100553");
  load(onia, kf_Upsilon_3S,"200553");
  load(onia, kf_chi_b0_1P, "10551");
  load(onia, kf_chi_b1_1P, "20553");
  load(onia, kf_chi_b2_1P, "555");
  load(onia, kf_chi_b0_2P, "110551");
  load(onia, kf_chi_b1_2P, "120553");
  load(onia, kf_chi_b2_2P, "100555");

  // --- Octets ---
  load(onia, kf_1S0_c_8_eta_c,        "1S0_441");
  load(onia, kf_1S0_c_8_J_psi_1S,     "1S0_443");
  load(onia, kf_1S0_c_8_psi_2S,       "1S0_100443");
  load(onia, kf_3S1_c_8_eta_c,        "3S1_441");
  load(onia, kf_3S1_c_8_J_psi_1S,     "3S1_443");
  load(onia, kf_3S1_c_8_psi_2S,       "3S1_100443");
  load(onia, kf_3S1_c_8_chi_c0_1P,    "3S1_10441");
  load(onia, kf_3S1_c_8_chi_c1_1P,    "3S1_20443");
  load(onia, kf_3S1_c_8_chi_c2_1P,    "3S1_445");
  load(onia, kf_3P0_c_8_J_psi_1S,     "3P0_443");
  load(onia, kf_3P0_c_8_psi_2S,       "3P0_100443");
  load(onia, kf_3P1_c_8_J_psi_1S,     "3P1_443");
  load(onia, kf_3P1_c_8_psi_2S,       "3P1_100443");
  load(onia, kf_3P2_c_8_J_psi_1S,     "3P2_443");
  load(onia, kf_3P2_c_8_psi_2S,       "3P2_100443");

  loaded = true;
  msg_Out() <<" LDME for Quarkonia has been loaded.\n";
}


double GetLDME(kf_code kfc)
{
  auto it = ldme_map.find(kfc);
  if (it == ldme_map.end()) {
    return 0.;
  }
  return it->second;
}

void WarnZeroLDME(kf_code kfc)
{
  msg_Error()
    << "Warning: LDME for kf_code "
    << kfc
    << " is zero. Skipping channel...\n";
}