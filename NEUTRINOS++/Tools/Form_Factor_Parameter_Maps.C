#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "ATOOLS/Phys/Flavour_Tags.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const cpl_info::code & cpl) {
  if (cpl==cpl_info::unknown)      s<<setw(12)<<"unknown";
  if (cpl==cpl_info::scalar)       s<<setw(12)<<"scalar";
  if (cpl==cpl_info::pseudoscalar) s<<setw(12)<<"pseudoscalar";
  if (cpl==cpl_info::vector)       s<<setw(12)<<"vector";
  if (cpl==cpl_info::axialvector)  s<<setw(12)<<"axialvector";
  if (cpl==cpl_info::tensor)       s<<setw(12)<<"tensor";
  return s;
}

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const ff_info & info) {
  s<<"     "<<info.m_cpl<<" ["<<info.m_type<<", "<<info.m_params.size()<<" parameters]: \n";
  s<<"                  ";
  if (info.m_type!=ff_type::none && info.m_type!=ff_type::unknown) {
    for (size_t i=0;i<info.m_params.size();i++) s<<info.m_params[i]<<" "; s<<"\n";
  }
  return s;
}
  
std::ostream & NEUTRINOS::operator<<(std::ostream & s, Form_Factor_Entry & entry) {
  for (map<kf_code, list<ff_info * > >::iterator eit=entry.GetAllEntries()->begin();
       eit!=entry.GetAllEntries()->end();eit++) { 
    s<<"   - "<<Flavour(eit->first).TexName()<<" mediated transitions:\n";
    for (list<ff_info * >::iterator ffit=eit->second.begin();
	 ffit!=eit->second.end();ffit++) {
      s<<(**ffit)<<"\n";
    }
  }
  return s;
}

Form_Factor_Entry::~Form_Factor_Entry() {
  while (!m_entries.empty()) {
    m_entries.begin()->second.erase(m_entries.begin()->second.begin(),
				    m_entries.begin()->second.end());
  }
}

void Form_Factor_Entry::Add(const kf_code & prop, list<ff_info * > & ffinfos) {
  if (m_entries.find(prop)!=m_entries.end()) {
    msg_Error()<<"Error in "<<METHOD<<": will override form factor information for "
	       <<prop<<" in the transition "<<m_pair.first<<" -> "<<m_pair.second<<"\n";
    m_entries[prop].erase(m_entries[prop].begin(),m_entries[prop].end());
  }
  m_entries[prop] = ffinfos;
}

Form_Factor_Parameter_Maps::Form_Factor_Parameter_Maps() {
  InitialiseMaps();
}

Form_Factor_Parameter_Maps::~Form_Factor_Parameter_Maps() {
  while (!empty()) {
    delete begin()->second;
    erase(begin());
  }      
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

void Form_Factor_Parameter_Maps::InitialiseMaps() {
  m_alltransitions = Settings::GetMainSettings()["TRANSITIONS"];
  InitialiseConstants();
  InitialiseFFMaps();
}

void Form_Factor_Parameter_Maps::InitialiseConstants() {
  m_parameters.clear();
  for (const auto& key: m_alltransitions["Constants"].GetKeys()) {
    m_parameters[key] = m_alltransitions["Constants"][key].SetDefault(-1.0).Get<double>();
  }
}

void Form_Factor_Parameter_Maps::InitialiseFFMaps() {
  for (auto& transition: m_alltransitions["Transitions"].GetKeys()) {
    InitialiseFormFactor(transition);
  }
}

bool Form_Factor_Parameter_Maps::InitialiseFormFactor(string & name) {
  vector<kf_code> kfcs;
  vector<double>  parameters;
  list<ff_info *> ffinfos;
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
      }
      for (auto & key: transition[cpl].GetKeys()) {
	if (key=="params") {
	  string params = transition[cpl][key].SetDefault("").Get<string>();
	  if (ExtractParameters(params,parameters)) {
	  }
	}
      }
      cpl_info::code cplinfo = CplInfo(cpl);
      ff_info * ffinfo = new ff_info(cplinfo,type,parameters.size());
      for (size_t i=0;i<parameters.size();i++) ffinfo->m_params[i] = parameters[i];
      ffinfos.push_back(ffinfo);
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

bool Form_Factor_Parameter_Maps::ExtractParameters(string & name,vector<double> & parameters) {
  string buffer = name;
  while (true) {
    while (buffer.length()>0 && (buffer[0]==' ' || buffer[0]=='\t')) buffer.erase(0,1);
    if (buffer.length()==0) return (parameters.size()>0);
    size_t pos(Min(buffer.find(' '),buffer.length()));
    parameters.push_back(ToType<double>(buffer.substr(0,pos)));
    buffer=buffer.substr(pos);
  }
}
