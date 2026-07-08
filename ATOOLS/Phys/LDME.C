#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/LDME.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <unordered_map>
#include <list> 


std::unordered_map<kf_code, double > ldme_sin_map;
std::map<kf_code, std::map<kf_code, double>> ldme_oct_map;
decay_channels d_1S0_c, d_3S1_c, d_3P0_c, d_3P1_c, d_3P2_c;
std::vector<decay_channels> d_alldecays_charm;

using namespace ATOOLS;

////////////////////////////////////////////////////////////////////
//   This is a helper file for Quarkonium production across sherpa. 
//   It gives the Long Distance Matrix Elements (LDME) and important
//   for the octet decays
////////////////////////////////////////////////////////////////////

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
void LoadLDMEMap(){
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
  d_alldecays_charm = {d_3S1_c, d_1S0_c, d_3P0_c, d_3P1_c, d_3P2_c};

  LoadLDME();
  auto onia = Settings::GetMainSettings()["QUARKONIA"];
  for(auto state=d_alldecays_charm.begin(); state!=d_alldecays_charm.end();state++){
    for(auto decay : state->decays) {
          std::string key = state->name + "_" + std::to_string(decay);
          ldme_oct_map[state->kfc][decay] = onia[key].Get<double>();
          msg_Debugging() << " ldme_oct_map[" <<state->kfc<<"]["<< decay <<"] = "<< ldme_oct_map[state->kfc][decay]<< '\n';
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
  onia["1S0_441"].SetDefault(0.); //still need tuning
  onia["1S0_443"].SetDefault(53.49E-03);
  onia["1S0_100443"].SetDefault(48.806E-03);
  onia["3S1_441"].SetDefault(0.);
  onia["3S1_443"].SetDefault(6.414E-03);
  onia["3S1_100443"].SetDefault(4.943E-03);
  onia["3S1_10441"].SetDefault(13.591E-03);
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
  double cf = 3.*3./(2.*M_PI); //we convert to LDMEs to keep it consistent
  onia["441"].SetDefault(cf*0.810);
  onia["443"].SetDefault(cf*0.810);
  onia["100443"].SetDefault(cf*0.529);
  onia["10441"].SetDefault(cf*0.075);
  onia["20443"].SetDefault(3.*cf*0.075);
  onia["445"].SetDefault(5.*cf*0.075);
  //BBBAR
  onia["551"].SetDefault(cf*6.477);
  onia["553"].SetDefault(cf*6.477);
  onia["100553"].SetDefault(3.*cf*.234);
  onia["200553"].SetDefault(3.*cf*2.474);
  onia["10551"].SetDefault(cf*1.417);
  onia["20553"].SetDefault(3.*cf*1.417);
  onia["555"].SetDefault(5.*cf*1.417);
  onia["110551"].SetDefault(cf*1.654);
  onia["120553"].SetDefault(3.*cf*1.654);
  onia["100555"].SetDefault(5.*cf*1.654);

  auto load = [&](auto& group, kf_code kfc, const char* key)
  {
    auto entry = group[key];
    ldme_sin_map[kfc] = entry.template Get<double>();
    if (entry.IsSetExplicitly())
      msgLDME(kfc, ldme_sin_map[kfc]);
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

const std::map<kf_code, double>& GetLDMEmap(kf_code kfc)
{
  LoadLDMEMap();

  static const std::map<kf_code, double> empty;
  auto it = ldme_oct_map.find(kfc);
  if (it == ldme_oct_map.end()){
    msg_Out()<<"No LDME map for "<<std::to_string(kfc); 
    return empty;
  }
  return it->second;
}

namespace{
  double GetSingletLDME(kf_code singlet_kfc)
  {
    LoadLDME();
    auto it = ldme_sin_map.find(singlet_kfc);
    if (it == ldme_sin_map.end()) {
      return 0.;
    }
    return it->second;
  }
}
////////////////////////////////////////////////////////////////////
//    The octet LDME is the sum of the LDMEs
//    of all of its "decay" channels
////////////////////////////////////////////////////////////////////

double GetLDME(kf_code kfc)
{
  double totalLDME = 0.;
  if(kfc/100000 == 99){
    LoadLDMEMap();
    const auto& decayLDMEmap = GetLDMEmap(kfc);
    for(const auto& decayprob : decayLDMEmap){
      totalLDME += decayprob.second;
    }
  }
  const double singletLDME = GetSingletLDME(kfc);
  if(IsZero(totalLDME)&&!IsZero(singletLDME)) totalLDME = singletLDME;
  else if(IsZero(totalLDME)&&IsZero(singletLDME)) WarnZeroLDME(kfc);
  return totalLDME;
}

////////////////////////////////////////////////////////////////////
//    We choose the decay channels relative to the LDMEs
////////////////////////////////////////////////////////////////////

kf_code GetQuarkoniumDecayChannel(kf_code octet){
  const auto& singlets = GetLDMEmap(octet);
  kf_code singlet_kfc;
  double totalprob = GetLDME(octet);
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