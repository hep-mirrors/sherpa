#include "MODEL/Main/Coupling_Data.H"

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;

void Coupling_Data::Calculate()
{
  if (p_scl==NULL) return;
  m_fac=(*p_cpl)(*p_scl)/m_def;
  msg_Debugging()<<METHOD<<"("<<this<<"): scl = "
		 <<sqrt(*p_scl)<<" -> "<<*this<<"\n";
}

void Coupling_Map::Calculate() const
{
  for (const_iterator cit(begin());
       cit!=end();++cit) cit->second->Calculate();
}

Coupling_Data *Coupling_Map::Get(const std::string &tag) const
{
  const_iterator cit(find(tag));
  if (cit!=end()) return cit->second;
  return NULL;
}

namespace MODEL {

  std::ostream &operator<<(std::ostream &str,const Coupling_Data &cd)
  {
    return str<<"'"<<cd.ID()<<"'{fac="
	      <<cd.Factor()<<",cpl="<<cd.Default()<<"}";
  }

}
