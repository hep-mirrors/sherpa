#include "COMIX/Phasespace/PS_Current.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

using namespace COMIX;
using namespace ATOOLS;

void PS_Current::SetSCC(Current_Base *const scc)
{
  p_scc=scc;
  m_psinfo="";
  m_psinfo=PSInfo();
  if (p_scc) m_psinfo+="_SC"+p_scc->PSInfo();
}

void PS_Current::ResetJ()
{
  m_j.resize(0); 
  m_zero=true;
}

void PS_Current::ResetZero()
{
  for (Vertex_Vector::const_iterator vit(m_out.begin());
       vit!=m_out.end();++vit) if (!(*vit)->Zero()) return;
  ResetJ();
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) (*vit)->SetZero();
}

void PS_Current::AddJ(const PS_Info &j)
{
  for (PS_Info_Vector::iterator cit(m_j.begin());
       cit!=m_j.end();++cit) if (j==*cit) return; 
  m_j.push_back(j);
  m_zero=false;
}

const std::vector<PS_Info> &PS_Current::J() const
{
  return m_j;
}

void PS_Current::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_id<<" {\n";
  msg_Indent();
#endif
  ResetJ();
  Vertex_Vector::const_iterator vit(m_in.begin());
  // calculate subcurrents
  for (;vit!=m_in.end();++vit) (*vit)->Evaluate();
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
  Print();
#endif
}

void PS_Current::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			    const int cr,const int ca)
{
  m_p=p;
  ResetJ();
  PS_Info i(cr,ca,1<<m_key);
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): '+' "<<m_id<<" "<<i<<" "<<m_fl<<"\n";
#endif
  AddJ(i);
}

void PS_Current::SetGauge(const ATOOLS::Vec4D &k)
{
}

void PS_Current::AddPropagator()
{
}

Complex PS_Current::Contract(const Current_Base &c,
			     const long unsigned int &hm,
			     const long unsigned int &hp) const
{
  return 0.0;
}

void PS_Current::Contract(const Current_Base &c,
			  const LongUIntMap_Matrix &hmp,
			  Complex_Vector &ress) const
{
}

char PS_Current::Type() const
{
  return 'P';
}   

void PS_Current::Print() const
{
  if (!msg_LevelIsDebugging()) return;
  std::string id(m_id.empty()?"<no entry>":ToString(m_id.front()));
  for (size_t i(1);i<m_id.size();++i) id+=","+ToString(m_id[i]);
  msg_Debugging()<<'['<<id<<"]{"<<m_id.size()<<","
		 <<m_key<<"}("<<m_oew<<","<<m_oqcd<<")("
		 <<(m_dir<0?Flav():Flav().Bar())<<")"
		 <<(m_dir==0?'-':m_dir>0?'I':'O')<<"{\n";
  if (m_p!=Vec4D()) msg_Debugging()<<"m_p  : "<<m_p<<"\n";
  if (!m_j.empty()) msg_Debugging()<<"m_j  :\n";
  {
    msg_Indent();
    for (size_t i(0);i<m_j.size();++i) 
      msg_Debugging()<<m_j[i]<<"\n";
  }
  if (!m_in.empty()) msg_Debugging()<<"m_in : ("<<m_in.size()<<")\n";
  {
    msg_Indent();
    for (size_t i(0);i<m_in.size();++i) 
      msg_Debugging()<<*m_in[i]<<"\n";
  }
  if (!m_out.empty()) msg_Debugging()<<"m_out: ("<<m_out.size()<<")\n";
  {
    msg_Indent();
    for (size_t i(0);i<m_out.size();++i) 
      msg_Debugging()<<*m_out[i]<<"\n";
  }
  msg_Debugging()<<"}\n";
}

namespace COMIX {

  template <> void Current_Base::AddJ(const PS_Info &j)
  { 
    return static_cast<PS_Current*>(this)->AddJ(j); 
  }
  
  template <> const std::vector<PS_Info> &Current_Base::J() const
  { 
    return static_cast<const PS_Current*>(this)->J(); 
  }

}
