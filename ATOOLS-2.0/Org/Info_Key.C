#include "Info_Key.H"

#include "Integration_Info.H"

using namespace ATOOLS;

Info_Key::Info_Key():
  p_info(NULL),
  m_valuekey(0),
  m_weightkey(0) {}

void Info_Key::Assign(const std::string name,const size_t doubles,
		      const size_t vectors,Integration_Info *const info)
{
  p_info=info;
  m_name=name;
  p_info->AssignKey(*this,doubles,vectors); 
}

Info_Key::~Info_Key()
{
  if (p_info!=NULL) p_info->ReleaseKey(*this);
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Info_Key &key)
{
  str<<"(\""<<key.m_name<<"\",\""<<key.m_info<<"\") -> "
     <<key.p_info->m_doubles[key.m_valuekey]<<" "
     <<key.p_info->m_vectors[key.m_valuekey]<<" => ("
     <<key.p_info->m_weights[key.m_valuekey][key.m_weightkey]<<")";
  return str;
}

