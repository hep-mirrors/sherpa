#include "ATOOLS/Org/Info_Key.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Smart_Pointer.C"

namespace ATOOLS { template class SP(Integration_Info); }

using namespace ATOOLS;

Integration_Info::Integration_Info() 
{
#ifdef USING__Threading
  pthread_mutex_init(&m_mtx,NULL);
#endif
}

Integration_Info::~Integration_Info() 
{
#ifdef USING__Threading
  pthread_mutex_destroy(&m_mtx);
#endif
  for (String_MapPair_Map::const_iterator kmit(m_keymap.begin());
       kmit!=m_keymap.end();++kmit) {
    const String_KeyPair_Map &keymap(kmit->second.second);
    for (String_KeyPair_Map::const_iterator kit=keymap.begin();
	 kit!=keymap.end();++kit) {
      const std::vector<Info_Key*> &keys(kit->second.second);
      for (std::vector<Info_Key*>::const_iterator it=keys.begin();
	   it!=keys.end();++it)
	(*it)->p_info=NULL;
    }
  }
}

void Integration_Info::ResetAll()
{
  for (size_t i(0);i<m_doubles.size();++i) {
    m_status[i]=si::reset;
#ifdef USING__safe_reset
    for (size_t j(0);j<m_doubles[i].size();++j) 
      m_doubles[i][j]=UNDEFINED_DOUBLE;
    for (size_t j(0);j<m_vectors[i].size();++j) 
      m_vectors[i][j]=UNDEFINED_VECTOR;
#endif
    for (size_t j(0);j<m_weights[i].size();++j) 
      m_weights[i][j]=UNDEFINED_WEIGHT;
  }
}

void Integration_Info::AssignKey(Info_Key &key,const size_t doubles,
				 const size_t vectors)
{
#ifdef USING__Threading
  pthread_mutex_lock(&m_mtx);
#endif
  if (key.p_info!=NULL) {
    if (key.p_info!=this) 
      msg_Error()<<METHOD<<"(): Key '"<<&key<<"' assigned to '"
		 <<key.p_info<<"'."<<std::endl;
    else
#ifdef DEBUG__Integration_Info
      msg_Debugging()<<METHOD<<"(): Key '"<<&key
		     <<"' already assigned."<<std::endl;
#endif
#ifdef USING__Threading
    pthread_mutex_unlock(&m_mtx);
#endif
    return;
  }
  if (m_keymap.find(key.m_name)==m_keymap.end()) {
    m_keymap[key.m_name]=SizeT_KeyMap_Pair(m_doubles.size(),
					   String_KeyPair_Map());
    m_doubles.push_back(Double_Container(doubles));
    m_vectors.push_back(Vector_Container(vectors));
    m_weights.push_back(Double_Container());
    m_status.push_back(si::idle);
  }
  key.m_valuekey=m_keymap[key.m_name].first;
  String_KeyPair_Map &keys(m_keymap[key.m_name].second);
  String_KeyPair_Map::iterator kit(keys.find(key.m_info));
  if (kit==keys.end()) {
    keys[key.m_info]=SizeT_KeyVector_Pair(m_weights[key.m_valuekey].size(),
					  Key_Vector());
    m_weights[key.m_valuekey].push_back(0.);
    kit=keys.find(key.m_info);
  }
  key.m_weightkey=kit->second.first;
  for (Key_Vector::iterator kvit(kit->second.second.begin());
       kvit!=kit->second.second.end();++kvit) if (*kvit==&key) {
#ifdef USING__Threading
      pthread_mutex_unlock(&m_mtx);
#endif
      return;
    }
  kit->second.second.push_back(&key);
  key.p_info=this;
#ifdef USING__Threading
  pthread_mutex_unlock(&m_mtx);
#endif
}

void Integration_Info::ReleaseKey(Info_Key &key)
{
#ifdef USING__Threading
  pthread_mutex_lock(&m_mtx);
#endif
  String_MapPair_Map::iterator vit(m_keymap.find(key.m_name));
  if (vit==m_keymap.end()) {
    msg_Error()<<METHOD<<"(): Name '"<<key.m_name
	       <<"' not found. Cannot release key "<<&key<<"."<<std::endl;
#ifdef USING__Threading
    pthread_mutex_unlock(&m_mtx);
#endif
    return;
  }
  String_KeyPair_Map &keys(vit->second.second);
  String_KeyPair_Map::iterator wit(keys.find(key.m_info));
  if (wit==keys.end()) {
    msg_Error()<<METHOD<<"(): Info '"<<key.m_info
	       <<"' not found. Cannot release key "<<&key<<"."<<std::endl;
#ifdef USING__Threading
    pthread_mutex_unlock(&m_mtx);
#endif
    return;
  }
  for (Key_Vector::iterator kit(wit->second.second.begin());
       kit!=wit->second.second.end();++kit) {
    if (*kit==&key) { 
      wit->second.second.erase(kit);
      key.p_info=NULL;
#ifdef USING__Threading
      pthread_mutex_unlock(&m_mtx);
#endif
      return;
    }
  }
  msg_Error()<<METHOD<<"(): Pointer '"<<&key
	     <<"' not found. Cannot release key."<<std::endl;
#ifdef USING__Threading
  pthread_mutex_unlock(&m_mtx);
#endif
}

std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Double_Container &doubles)
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
 
std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Vector_Container &vectors)
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

std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Integration_Info &info)
{
  str<<"Integration_Info("<<&info<<") {\n";
  for (size_t i=0;i<info.m_doubles.size();++i) {
    str<<"  (*this)["<<i<<"] = "<<info.m_doubles[i]<<" "
       <<info.m_vectors[i]<<" => "<<info.m_weights[i]
       <<" => ("<<info.m_status[i]<<")\n";
  }
  str<<"}";
  return str;
}
