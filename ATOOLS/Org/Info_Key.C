#include "ATOOLS/Org/Info_Key.H"

#include "ATOOLS/Org/Integration_Info.H"

using namespace ATOOLS;

Info_Key::Info_Key():
  p_info(NULL),
  m_valuekey(0),
  m_weightkey(0) {}

void Info_Key::Assign(const std::string name,const size_t doubles,
		      const size_t vectors,Integration_Info *const info)
{
  m_name=name;
  info->AssignKey(*this,doubles,vectors);
}

void Info_Key::Assign(const std::string name,const size_t doubles,
		      const size_t vectors,const SP(Integration_Info) &info)
{
  m_name=name;
  info->AssignKey(*this,doubles,vectors);
}

Info_Key::~Info_Key()
{
  if (p_info!=NULL) p_info->ReleaseKey(*this);
}

void Info_Key::SetInfo(const std::string info)
{
  Integration_Info *cinfo(p_info);
  if (cinfo!=NULL) cinfo->ReleaseKey(*this);
  m_info=info;
  if (cinfo!=NULL) cinfo->AssignKey(*this,0,0);
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Info_Key &key)
{
  str<<"(\""<<key.m_name<<"\",\""<<key.m_info<<"\") -> ";
  if (key.p_info==NULL) return str<<"NULL";
  return str<<key.p_info->m_doubles[key.m_valuekey]<<" "
	    <<key.p_info->m_vectors[key.m_valuekey]<<" => ("
	    <<key.p_info->m_weights[key.m_valuekey][key.m_weightkey]<<")";
}

