#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Remnants_Parameters* REMNANTS::rempars = NULL;

void pkparams::Output() {
  msg_Debugging()<<"-------------------------------------------------\n"
	   <<"Form = "<<m_form<<", Recoil = "<<m_recoil<<"\n";
  if (m_form!=pkform::none) {
    for (map<string, double>::iterator pit=m_params.begin();pit!=m_params.end();pit++) {
      msg_Debugging()<<std::left<<std::setw(32)<<pit->first<<": "<<pit->second<<"\n";
    }
  }
  msg_Debugging()<<"-------------------------------------------------\n";
}


Remnants_Parameters::Remnants_Parameters() {
  Initialize();
}

Remnants_Parameters::~Remnants_Parameters() {
  for (map<kf_code,pkparams *>::iterator it=m_params.begin();it!=m_params.end();it++) {
    delete it->second;
  }
  m_params.clear();
}

void Remnants_Parameters::Initialize() {
  pkparams * proton = new pkparams;
  m_params[2212]    = proton;
  proton->m_form    = pkform::gauss_limited;
  proton->m_recoil  = pkrecoil::beam_vs_shower; //pkrecoil::democratic;
  proton->m_params["SHOWER_INITIATOR_MEAN"]   = 1.0;
  proton->m_params["SHOWER_INITIATOR_SIGMA"]  = 1.1;
  proton->m_params["SHOWER_INITIATOR_Q2"]     = 0.77;
  proton->m_params["SHOWER_INITIATOR_KTMAX"]  = 2.7;
  proton->m_params["SHOWER_INITIATOR_KTEXPO"] = 5.12;
  proton->m_params["REFERENCE_ENERGY"]        = 7000.0;
  proton->m_params["ENERGY_SCALING_EXPO"]     = 0.08;
  proton->m_params["BEAM_SPECTATOR_MEAN"]     = 0.0;
  proton->m_params["BEAM_SPECTATOR_SIGMA"]    = 0.25;
  proton->m_params["BEAM_SPECTATOR_Q2"]       = 0.77;
  proton->m_params["BEAM_SPECTATOR_KTMAX"]    = 1.0;
  proton->m_params["BEAM_SPECTATOR_KTEXPO"]   = 5.0;
  // As default settings copy the proton defaults
  m_params[0]        = new pkparams(*proton);
}

void Remnants_Parameters::ReadParameters() {
  auto kperplist = Settings::GetMainSettings()["INTRINSIC_KPERP"];
  for (const auto& kfit : kperplist.GetKeys()) {
    pkparams * params = NULL;
    kf_code kfcode = ToType<kf_code>(kfit);
    if (s_kftable.find(kfcode)==s_kftable.end()) continue;
    if (m_params.find(kfcode)==m_params.end()) {
      msg_Debugging()<<"* Reading intrinsic kperp parameters for new particle: "<<kfcode<<"\n";
      params           = new pkparams;
      m_params[kfcode] = params;
    }
    else {
      msg_Debugging()<<"* Updating intrinsic kperp parameters for: "<<kfcode<<"\n";
      params = m_params.find(kfcode)->second;
    }
    for (const auto& pname : kperplist[kfit].GetKeys()) {
      auto s = kperplist[kfit][pname];
      if (pname=="FORM")
	params->m_form = SelectForm(s.SetDefault("none").Get<std::string>());
      else if (pname=="RECOIL")
	params->m_recoil = SelectRecoil(s.SetDefault("beam_vs_shower").Get<std::string>());
      else
	params->m_params[pname] = s.SetDefault(-1.).Get<double>();
    }
  }
  msg_Debugging()<<"============================================================\n"
	   <<METHOD<<" results in final parameters:\n";
  for (map<long unsigned int, pkparams * >::iterator pit=m_params.begin();
       pit!=m_params.end();pit++) {
    msg_Debugging()<<"* "<<pit->first<<"\n";
    pit->second->Output();
  }
  msg_Debugging()<<"============================================================\n";
}

pkform::code Remnants_Parameters::SelectForm(const std::string & form) {
  pkform::code pkf = pkform::undefined;
  if (form=="None" || form=="none") pkf = pkform::none;
  else if (form=="gauss")           pkf = pkform::gauss;
  else if (form=="gauss_limited")   pkf = pkform::gauss_limited;
  else if (form=="dipole")          pkf = pkform::dipole;
  else if (form=="dipole_limited")  pkf = pkform::dipole_limited;
  else THROW(not_implemented,"Intrinsic KPerp model form not implemented.");
  return pkf;
}

pkrecoil::code Remnants_Parameters::SelectRecoil(const std::string & form) {
  pkrecoil::code pkr = pkrecoil::undefined;
  if (form=="democratic")          pkr = pkrecoil::democratic;
  else if (form=="beam_vs_shower") pkr = pkrecoil::beam_vs_shower;
  else THROW(not_implemented,"Intrinsic KPerp recoil model not implemented.")
  return pkr;
}


const pkform::code & Remnants_Parameters::GetForm(const Flavour & beamflav) {
  map<long unsigned int,pkparams *>::iterator fit = m_params.find(beamflav.Kfcode());
  return (fit!=m_params.end() ? fit->second->m_form : m_params[0]->m_form);
}

const pkrecoil::code & Remnants_Parameters::GetRecoil(const Flavour & beamflav) {
  map<long unsigned int,pkparams *>::iterator fit = m_params.find(beamflav.Kfcode());
  return (fit!=m_params.end() ? fit->second->m_recoil : m_params[0]->m_recoil);
}

const double & Remnants_Parameters::operator()(const Flavour & beamflav,const string & tag) {
  map<long unsigned int,pkparams *>::iterator fit = m_params.find(beamflav.Kfcode());
  pkparams * pkparams = (fit!=m_params.end() ? fit->second : m_params[0]);
  map<string, double>::iterator tit = pkparams->m_params.find(tag);
  if (tit==pkparams->m_params.end()) THROW(fatal_error,tag+" not found in remnants parameters.")
  return tit->second;
}

