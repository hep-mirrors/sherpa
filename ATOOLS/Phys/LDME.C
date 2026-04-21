#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/LDME.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <unordered_map>
#include <list> 


std::unordered_map<kf_code, double > ldme_map;
std::map<kf_code, std::map<kf_code, double>> LDMEmap;
decay_channels d_1S0_c, d_3S1_c, d_3P0_c, d_3P1_c, d_3P2_c;
std::list<decay_channels> d_ch_c;
std::map<kf_code, double> singlets;

using namespace ATOOLS;

void WarnZeroLDME(kf_code kfc)
{
  msg_Error()
    << "Warning: LDME for kf_code "
    << kfc
    << " is zero. Skipping channel...\n";
}

void msgLDME(kf_code kfc, double LDME)
{
  msg_Info()
    <<"LDME for kf_code "
    << kfc
    << " has been set to "
    << LDME<<"\n";
}
void LoadLDMEmap(){
  static bool loaded = false;
  if(loaded) return;

  d_1S0_c.kfc = kf_1S0_c; 
  d_1S0_c.name = "1S0"; 
  d_1S0_c.decays = {kf_J_psi_1S, kf_eta_c_1S, kf_psi_2S};
  d_3S1_c.kfc = kf_3S1_c; 
  d_3S1_c.name = "3S1"; 
  d_3S1_c.decays = {kf_J_psi_1S, kf_eta_c_1S, kf_psi_2S, kf_chi_c0_1P, kf_chi_c1_1P, kf_chi_c2_1P};
  d_3P0_c.kfc = kf_3P0_c; 
  d_3P0_c.name = "3P0";
  d_3P0_c.decays = {kf_J_psi_1S, kf_psi_2S};
  d_3P1_c.kfc = kf_3P1_c; 
  d_3P1_c.name = "3P1";
  d_3P1_c.decays = {kf_J_psi_1S, kf_psi_2S};
  d_3P2_c.kfc = kf_3P2_c; 
  d_3P2_c.name = "3P2";
  d_3P2_c.decays = {kf_J_psi_1S, kf_psi_2S};
  d_ch_c = {d_3S1_c, d_1S0_c, d_3P0_c, d_3P1_c, d_3P2_c};

  LoadLDME();
  auto onia = Settings::GetMainSettings()["QUARKONIA"];
  for(auto state=d_ch_c.begin(); state!=d_ch_c.end();state++){
    for(auto decay : state->decays) {
          std::string key = state->name + "_" + std::to_string(decay);
          LDMEmap[state->kfc][decay] = onia[key].Get<double>();
          msg_Debugging() << " LDMEmap[" <<state->kfc<<"]["<< decay <<"] = "<< LDMEmap[state->kfc][decay]<< '\n';
    }
  }
  loaded = true;
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

  loaded = true;
  msg_Debugging() <<" LDME for Quarkonia has been loaded.\n";
}


double GetLDME(kf_code kfc)
{
  LoadLDME();
  auto it = ldme_map.find(kfc);
  if (it == ldme_map.end()) {
    return 0.;
  }
  return it->second;
}

std::map<kf_code, double> GetLDMEmap(kf_code kfc)
{
  LoadLDMEmap();
  return LDMEmap[kfc];
}

double GetTotalLDME(kf_code kfc)
{
  LoadLDMEmap();
  std::map<kf_code, double> decayLDMEmap = GetLDMEmap(kfc);
  double totalLDME = 0.;
  for(auto decayprob : decayLDMEmap){
    totalLDME += decayprob.second;
  }
  if(IsZero(totalLDME)&&!IsZero(GetLDME(kfc))) totalLDME = GetLDME(kfc);
  else if(IsZero(totalLDME)&&IsZero(GetLDME(kfc))) WarnZeroLDME(kfc);
  return totalLDME;
}

kf_code GetChannel(kf_code octet){
  singlets = GetLDMEmap(octet);
  kf_code singlet_kfc;
  double totalprob = GetTotalLDME(octet);
  if(IsZero(totalprob)) return octet;
  double r = ran->Get();
  double decayProb = 0.;
  for (auto it = singlets.begin(); it !=singlets.end(); it++){
    decayProb += it->second/totalprob;
      if(r <= decayProb) {
        singlet_kfc = it->first;
        return singlet_kfc;
      }
  }
  msg_Error()<< METHOD <<" Decay channel for octet "<< octet<<" not found.\n";
  return octet;
}