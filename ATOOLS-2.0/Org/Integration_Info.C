#include "Info_Key.H"

#include "Message.H"

using namespace ATOOLS;

Integration_Info::Integration_Info() 
{
}

void Integration_Info::ResetAll()
{
  for (size_t i=0;i<m_doubles.size();++i) {
    m_status[i]=si::reset;
#ifdef USING__safe_reset
    for (size_t j=0;j<m_doubles[i].size();++j) m_doubles[i][j]=UNDEFINED_DOUBLE;
    for (size_t j=0;j<m_vectors[i].size();++j) m_vectors[i][j]=UNDEFINED_VECTOR;
#endif
    for (size_t j=0;j<m_weights[i].size();++j) m_weights[i][j]=UNDEFINED_WEIGHT;
  }
}

void Integration_Info::AssignKey(Info_Key &key,const size_t doubles,const size_t vectors)
{
  if (m_keymap.find(key.m_name)==m_keymap.end()) {
    m_keymap[key.m_name]=std::pair<size_t,std::map<std::string,
      std::pair<size_t,std::vector<Info_Key*> > > >
      (m_doubles.size(),std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > >());
    m_doubles.push_back(Double_Container(doubles));
    m_vectors.push_back(Vector_Container(vectors));
    m_weights.push_back(Double_Container());
    m_status.push_back(si::idle);
  }
  key.m_valuekey=m_keymap[key.m_name].first;
  std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > > &keys=m_keymap[key.m_name].second;
  std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > >::iterator kit=keys.find(key.m_info);
  if (kit==keys.end()) {
    keys[key.m_info]=std::pair<size_t,std::vector<Info_Key*> >(m_weights[key.m_valuekey].size(),
							       std::vector<Info_Key*>());
    m_weights[key.m_valuekey].push_back(0.);
    kit=keys.find(key.m_info);
  }
  key.m_weightkey=kit->second.first;
  kit->second.second.push_back(&key);
}

void Integration_Info::ReleaseKey(Info_Key &key)
{
  std::map<std::string,std::pair<size_t,
    std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > > > >::iterator 
    vit=m_keymap.find(key.m_name);
  if (vit==m_keymap.end()) return;
  size_t vpos=vit->second.first;
  std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > > &keys=vit->second.second;
  std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > >::iterator 
    wit=keys.find(key.m_info);
  if (wit==keys.end()) return;
  size_t wpos=wit->second.first;
  for (std::vector<Info_Key*>::iterator kit=wit->second.second.begin();
       kit!=wit->second.second.end();++kit) {
    if (*kit==&key) { 
      wit->second.second.erase(kit);
      return;
    }
  }
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Double_Container &doubles)
{
  if (doubles.size()==0) return str<<"{<no entries>}";
#ifdef __GNUC__
#if __GNUC__ > 2
  std::ios_base::fmtflags flags=str.flags();
#else
  std::ios::fmtflags flags=str.flags();
#endif
#else
  std::ios::fmtflags flags=str.flags();
#endif
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
 
std::ostream &ATOOLS::operator<<(std::ostream &str,const Vector_Container &vectors)
{
  if (vectors.size()==0) return str<<"{<no entries>}";
#ifdef __GNUC__
#if __GNUC__ > 2
  std::ios_base::fmtflags flags=str.flags();
#else
  std::ios::fmtflags flags=str.flags();
#endif
#else
  std::ios::fmtflags flags=str.flags();
#endif
  str.precision(6);
  str<<"{";
  for (size_t i=0;i<vectors.size();++i) {
    str<<vectors[i]<<",";
  }
  str<<"\b}";
  str.setf(flags);
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Integration_Info &info)
{
  str<<"Integration_Info("<<&info<<") {\n";
  for (size_t i=0;i<info.m_doubles.size();++i) {
    str<<"  (*this)["<<i<<"] = "<<info.m_doubles[i]<<" "<<info.m_vectors[i]<<" => "
       <<info.m_weights[i]<<" => ("<<info.m_status[i]<<")\n";
  }
  str<<"}";
  return str;
}
