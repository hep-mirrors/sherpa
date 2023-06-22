#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "ATOOLS/Phys/Flavour_Tags.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

namespace NEUTRINOS {
  Form_Factor_Parameter_Maps * ffs(NULL);
}

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

std::ostream & NEUTRINOS::operator<<(std::ostream & s, Form_Factor_Entry & entry) {
  for (map<kf_code, map<cpl_info::code, ff_info * > >::iterator eit=entry.GetAll()->begin();
       eit!=entry.GetAll()->end();eit++) { 
    s<<"   - "<<Flavour(eit->first).TexName()<<" mediated transitions:\n";
    for (map<cpl_info::code, ff_info * > ::iterator ffit=eit->second.begin();
	 ffit!=eit->second.end();ffit++) {
      s<<(*ffit->second)<<"\n";
    }
  }
  return s;
}

Form_Factor_Entry::~Form_Factor_Entry() {
  while (!m_entries.empty()) {
    map<cpl_info::code, ff_info * > ffinfos = m_entries.begin()->second;  
    while (!ffinfos.empty()) {
      delete (ffinfos.begin()->second); ffinfos.erase(ffinfos.begin());
    }
    m_entries.erase(m_entries.begin());
  }
}

void Form_Factor_Entry::Add(const kf_code & prop, map<cpl_info::code, ff_info * > & ffinfos) {
  if (m_entries.find(prop)!=m_entries.end()) {
    msg_Error()<<"Error in "<<METHOD<<": will override form factor information for "
	       <<prop<<" in the transition "<<m_pair.first<<" -> "<<m_pair.second<<"\n";
    map<cpl_info::code, ff_info * > ffinfos = m_entries[prop];
    while (!ffinfos.empty()) {
      delete (ffinfos.begin()->second); ffinfos.erase(ffinfos.begin());
    }
  }
  m_entries[prop] = ffinfos;
}

Form_Factor_Parameter_Maps::Form_Factor_Parameter_Maps() {}

Form_Factor_Parameter_Maps::~Form_Factor_Parameter_Maps() {
  while (!empty()) {
    delete begin()->second;
    erase(begin());
  }      
}

Form_Factor_Base *
Form_Factor_Parameter_Maps::GetFF(kf_code & n1,kf_code & n2,kf_code & prop,cpl_info::code & cpl) {
  ff_info * info1 = FindInfo(n1,n2,prop,cpl);
  ff_info * info2 = FindInfo(n2,n1,prop,cpl);

  ff_info * info;
  if (n1 != n2) {
    if (info1 == NULL && info2 != NULL) info = info2;
    else if (info2 == NULL && info1 != NULL) info = info1;
    else info = NULL;
  } else info = info1;

  if (info!=NULL) {
    switch (info->m_type) {
    case ff_type::constant:
      return new Constant_Form_Factor(*info);
    case ff_type::dipole:
      return new Dipole_Form_Factor(*info);
    case ff_type::lambda_dipole:
      return new Lambda_Dipole_Form_Factor(*info);
    case ff_type::neutron_electric:
      return new Neutron_Electric_Form_Factor(*info);
    case ff_type::exponential:
      return new Exponential_Form_Factor(*info);
    case ff_type::Gaussian:
      return new Gaussian_Form_Factor(*info);
    case ff_type::Kelly:
      return new Kelly_Form_Factor(*info);
    case ff_type::BBBA:
      return new BBBA_Form_Factor(*info);
    case ff_type::ArringtonHill:
      return new ArringtonHill_Form_Factor(*info);
    case ff_type::Helm:
      return new Helm_Form_Factor(*info);
    case ff_type::Lovato:
      return new Lovato_Form_Factor(*info);
    case ff_type::unknown:
      msg_Error()<<"Error in "<<METHOD<<": found unknown FF type.\n"
		 <<"    Will return Dummy Form Factor for "
		 <<Flavour(n1)<<" -> "<<Flavour(n2)<<" + "<<Flavour(prop)<<"\n";
    case ff_type::none:
    default:
      break;
    }
  }
  //return new Dummy_Form_Factor(cpl);
  // JW: Does it make more sense to return a 0 form factor for the form factor type if not defined in runcard?
  return new Zero_Form_Factor(cpl);
}

double Form_Factor_Parameter_Maps::GetModelParms(std::string constantName) {
    for (map<string, double>::iterator pit=m_parameters.begin();pit!=m_parameters.end();pit++) {
      if (pit->first.c_str() == constantName) return pit->second;
    }
    //Return non-valid ckm transition
    if ( constantName == "Vuu_Vdd") return 0.;
    if ( constantName == "Vqq") return 1.;

    return 0.0;
}

bool Form_Factor_Parameter_Maps::ContainsFormFactorType(cpl_info::code & cpl) {
  kf_code prop_kf_P   = kf_photon;
  kf_code prop_kf_Z   = kf_Z;
  kf_code prop_kf_W   = kf_Wplus;

  for (std::map<std::pair<kf_code, kf_code > , Form_Factor_Entry * >::iterator it=begin(); it!=end();it++) {
      kf_code N1 = it->first.first;
      kf_code N2 = it->first.second;

      ff_info * is_QED = FindInfo(N1,N2,prop_kf_P,cpl);
      ff_info * is_Weak_NC = FindInfo(N1,N2,prop_kf_Z,cpl);
      ff_info * is_Weak_CC = FindInfo(N1,N2,prop_kf_W,cpl);
      if (is_QED != NULL || is_Weak_NC != NULL || is_Weak_CC != NULL) return true;
  }
  return false;
}

bool Form_Factor_Parameter_Maps::ContainsFormFactorTransition(kf_code & pid1, kf_code & pid2) {
  for (std::map<std::pair<kf_code, kf_code > , Form_Factor_Entry * >::iterator it=begin(); it!=end();it++) {
      kf_code N1 = it->first.first;
      kf_code N2 = it->first.second;
      if ( (N1 != pid1 & N2 != pid2) || (N1 != pid2 & N2 != pid1)) return true;
  }
  return false;
}

void Form_Factor_Parameter_Maps::Output() {
  msg_Info()<<"--------------------------------------------------------------------------------------\n";
  msg_Info()<<"* Global constants and parameters:\n";
  for (map<string, double>::iterator pit=m_parameters.begin();pit!=m_parameters.end();pit++) {
    msg_Info()<<"   "<<setw(12)<<pit->first.c_str()<<": "<<pit->second<<"\n";
  }
  for (std::map<std::pair<kf_code, kf_code > , Form_Factor_Entry * >::iterator it=begin();
       it!=end();it++) {
    msg_Info()<<"* Transitions for "<<Flavour(it->first.first)<<" --> "<<Flavour(it->first.second)<<":\n"
	      <<(*it->second);;
  }
  msg_Info()<<"--------------------------------------------------------------------------------------\n";
}

void Form_Factor_Parameter_Maps::Initialize() {
  m_alltransitions = Settings::GetMainSettings()["TRANSITIONS"];
  InitializeConstants();
  InitializeFFMaps();
}

void Form_Factor_Parameter_Maps::InitializeConstants() {
  m_parameters.clear();
  for (const auto& key: m_alltransitions["Constants"].GetKeys()) {
    m_parameters[key] = m_alltransitions["Constants"][key].SetDefault(-1.0).Get<double>();
  }
}

void Form_Factor_Parameter_Maps::InitializeFFMaps() {
  for (auto& transition: m_alltransitions["Transitions"].GetKeys()) {
    InitializeFormFactor(transition);
  }
}

bool Form_Factor_Parameter_Maps::InitializeFormFactor(string & name) {
  vector<kf_code> kfcs;
  vector<double>  parameters;
  map<cpl_info::code, ff_info *> ffinfos;
  Scoped_Settings transition = m_alltransitions["Transitions"][name]; 
  if (ExtractFlavours(name,kfcs)) {
    pair<kf_code, kf_code> kfpair = make_pair(kfcs[0], kfcs[1]);
    if (find(kfpair)==end()) (*this)[kfpair] = new Form_Factor_Entry(kfpair);
    for (auto & cpl: transition.GetKeys()) {
      parameters.clear();
      ff_type::code type = ff_type::none;
      for (auto & key: transition[cpl].GetKeys()) {
        if (key=="form") {
          string form = transition[cpl][key].SetDefault("none").Get<string>();
          type        = FFType(form);
        }
        parameters  = transition[cpl]["params"].SetDefault(vector<double>{}).GetVector<double>();
      }
      cpl_info::code cplinfo = CplInfo(cpl);
      ff_info * ffinfo = new ff_info(cplinfo,type,parameters.size());
      for (size_t i=0;i<parameters.size();i++) ffinfo->m_params[i] = parameters[i];
      ffinfos[cplinfo] = ffinfo;
    }
    Form_Factor_Entry * entry = (*this)[kfpair];
    entry->Add(kfcs[2], ffinfos);
    return true;
  }
  return false;
}


bool Form_Factor_Parameter_Maps::ExtractFlavours(string & name,vector<kf_code> & kfcs) {
  string buffer = name;
  while (true) {
    while (buffer.length()>0 && (buffer[0]==' ' || buffer[0]=='\t')) buffer.erase(0,1);
    if (buffer.length()==0) return (kfcs.size()>0);
    size_t pos(Min(buffer.find(' '),buffer.length()));
    kfcs.push_back(ToType<kf_code>(buffer.substr(0,pos)));
    buffer=buffer.substr(pos);
  }
}


