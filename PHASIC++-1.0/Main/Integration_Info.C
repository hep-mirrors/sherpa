#include "Integration_Info.H"

#include "Message.H"

using namespace PHASIC;

Info_Key::Info_Key():
  p_info(NULL),
  m_valuekey(0),
  m_weightkey(0) {}

void Info_Key::Assign(const std::string name,const size_t doubles,
		      const size_t vectors,Integration_Info *info)
{
  p_info=info;
  m_name=name;
  p_info->AssignKey(*this,doubles,vectors); 
}

Info_Key::~Info_Key()
{
  p_info->ReleaseKey(*this);
}

Integration_Info::Integration_Info() 
{
}

void Integration_Info::ResetAll()
{
  for (size_t i=0;i<m_doubles.size();++i) {
    m_status[i]=si::reset;
    // for (size_t j=0;j<m_doubles[i].size();++j) m_doubles[i][j]=UNDEFINED_DOUBLE;
    // for (size_t j=0;j<m_vectors[i].size();++j) m_vectors[i][j]=UNDEFINED_VECTOR;
    for (size_t j=0;j<m_weights[i].size();++j) m_weights[i][j]=UNDEFINED_WEIGHT;
  }
}

void Integration_Info::AssignKey(Info_Key &key,const size_t doubles,const size_t vectors)
{
  if (m_keymap.find(key.m_name)==m_keymap.end()) {
    m_keymap[key.m_name]=std::pair<size_t,std::vector<Info_Key*> >
      (m_doubles.size(),std::vector<Info_Key*>());
    m_doubles.push_back(Double_Container(doubles));
    m_vectors.push_back(Vector_Container(vectors));
    m_weights.push_back(Double_Container());
    m_status.push_back(si::idle);
  }
  key.m_valuekey=m_keymap[key.m_name].first;
  std::vector<Info_Key*> &keys=m_keymap[key.m_name].second;
  size_t i;
  for (i=0;i<keys.size();++i) {
    if (keys[i]->m_info==key.m_info) { 
      key.m_weightkey=keys[i]->m_weightkey;
      break;
    }
  }
  if (i==keys.size()) {
    key.m_weightkey=keys.size();
    keys.push_back(&key);
    m_weights[key.m_valuekey].push_back(0.);
  }
  ATOOLS::msg.Tracking()<<"key mapping: ("<<key.m_valuekey<<","<<key.m_weightkey<<") <=> (\""
			<<key.m_name<<"\",\""<<key.m_info<<"\")\n";
}

void Integration_Info::ReleaseKey(Info_Key &key)
{
}

std::ostream &PHASIC::operator<<(std::ostream &str,const Double_Container &doubles)
{
  if (doubles.size()==0) return str<<"{<no entries>}";
  std::ios_base::fmtflags flags=str.flags();
  str.precision(6);
  str<<"{";
  for (size_t i=0;i<doubles.size();++i) {
    // str.width(13); 
    str<<doubles[i]<<",";
  }
  str<<"\b}";
  str.setf(flags);
  return str;
}
 
std::ostream &PHASIC::operator<<(std::ostream &str,const Vector_Container &vectors)
{
  if (vectors.size()==0) return str<<"{<no entries>}";
  std::ios_base::fmtflags flags=str.flags();
  str.precision(6);
  str<<"{";
  for (size_t i=0;i<vectors.size();++i) {
    str<<vectors[i]<<",";
  }
  str<<"\b}";
  str.setf(flags);
  return str;
}

std::ostream &PHASIC::operator<<(std::ostream &str,const Info_Key &key)
{
  str<<"(\""<<key.m_name<<"\",\""<<key.m_info<<"\") -> "
     <<key.p_info->m_doubles[key.m_valuekey]<<" => ("
     <<key.p_info->m_weights[key.m_valuekey][key.m_weightkey]<<")";
  return str;
}

std::ostream &PHASIC::operator<<(std::ostream &str,const Integration_Info &info)
{
  str<<"Integration_Info("<<&info<<") {\n";
  for (size_t i=0;i<info.m_doubles.size();++i) {
    str<<"  (*this)["<<i<<"] = "<<info.m_doubles[i]<<" "<<info.m_vectors[i]<<" => "
       <<info.m_weights[i]<<" => ("<<info.m_status[i]<<")\n";
  }
  str<<"}";
  return str;
}
