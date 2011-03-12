#include "COMIX/Currents/Current.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

using namespace COMIX;
using namespace ATOOLS;

template <typename CType>
void Current<CType>::ResetJ()
{
  m_j.resize(0);
  m_zero=true;
}

template <typename CType>
void Current<CType>::ResetZero()
{
  for (Vertex_Vector::const_iterator vit(m_out.begin());
       vit!=m_out.end();++vit) if (!(*vit)->Zero()) return;
  ResetJ();
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) (*vit)->SetZero();
}

template <typename CType>
void Current<CType>::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_id<<" {\n";
  msg_Indent();
#endif
  ResetJ();
  Vertex_Vector::const_iterator vit(m_in.begin());
  // calculate outgoing momentum
  m_p=(*vit)->JA()->P()+(*vit)->JB()->P();
  // calculate subcurrents
  for (;vit!=m_in.end();++vit) (*vit)->Evaluate();
  if (!m_out.empty() && !m_zero) AddPropagator();
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
  Print();
#endif
}

template <typename CType>
void Current<CType>::Print() const
{
  if (!msg_LevelIsDebugging()) return;
  std::string id(m_id.empty()?"<no entry>":ToString(m_id.front()));
  for (size_t i(1);i<m_id.size();++i) id+=","+ToString(m_id[i]);
  msg_Debugging()<<'['<<id<<"]"<<m_fid<<"{"<<m_id.size()<<","
		 <<m_key<<"}("<<m_oew<<","<<m_oqcd<<"|"<<m_ntc<<")("
		 <<(m_dir<=0?Flav():Flav().Bar())<<")"
		 <<(m_dir==0?"":m_dir>0?"I":"O")
		 <<(m_cut?"C"+ToString(m_cut):"")<<"{\n";
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

template <typename CType> void Current_Base::AddJ(const CType &j)
{ 
  typename std::vector<CType> 
    &cj(static_cast<Current<CType>*>(this)->m_j); 
  for (typename std::vector<CType>::iterator cit(cj.begin());
       cit!=cj.end();++cit)
    if (j.H(0)==cit->H(0) && j.H(1)==cit->H(1) &&
	j(0)==(*cit)(0) && j(1)==(*cit)(1)) { 
      *cit+=j;
      return; 
    }
  cj.push_back(j);
  m_zero=false;
}

template <typename CType> const std::vector<CType> &Current_Base::J() const
{ 
  return static_cast<const Current<CType>*>(this)->m_j; 
}

namespace COMIX {

  template void Current_Base::AddJ(const DCScalar &j);
  template const std::vector<DCScalar> &Current_Base::J() const;

  template class Current<DCScalar>;
  
  template void Current_Base::AddJ(const DDSpinor &j);
  template const std::vector<DDSpinor> &Current_Base::J() const;
  
  template class Current<DDSpinor>;
  
  template void Current_Base::AddJ(const DCVec4D &j);
  template const std::vector<DCVec4D> &Current_Base::J() const;
  
  template class Current<DCVec4D>;
  
  template void Current_Base::AddJ(const DCAsT4D &j);
  template const std::vector<DCAsT4D> &Current_Base::J() const;
  
  template class Current<DCAsT4D>;


  template void Current_Base::AddJ(const QCScalar &j);
  template const std::vector<QCScalar> &Current_Base::J() const;

  template class Current<QCScalar>;

  template void Current_Base::AddJ(const QDSpinor &j);
  template const std::vector<QDSpinor> &Current_Base::J() const;

  template class Current<QDSpinor>;

  template void Current_Base::AddJ(const QCVec4D &j);
  template const std::vector<QCVec4D> &Current_Base::J() const;

  template class Current<QCVec4D>;

  template void Current_Base::AddJ(const QCAsT4D &j);
  template const std::vector<QCAsT4D> &Current_Base::J() const;

  template class Current<QCAsT4D>;

}
