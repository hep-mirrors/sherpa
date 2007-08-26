#include "Current.H"

#include "Message.H"
#include "MyStrStream.H"
#include "STL_Tools.H"

using namespace EXTRAXS;
using namespace ATOOLS;

template <typename CType>
void Current<CType>::ResetJ()
{
  m_j.resize(0); 
  for (Vertex_Vector::const_iterator vit=m_out.begin();
       vit!=m_out.end();++vit) (*vit)->SetZero(false);
}

template <typename CType>
void Current<CType>::AddJ(const CType &j)
{
  for (typename std::vector<CType>::iterator cit(m_j.begin());
       cit!=m_j.end();++cit)
    if (j(0)==(*cit)(0) && j(1)==(*cit)(1) && 
	j.H(0)==cit->H(0) && j.H(1)==cit->H(1)) { 
      *cit+=j;
      return; 
    }
  m_j.push_back(j);
}

template <typename CType>
const std::vector<CType> &Current<CType>::J() const
{
  return m_j;
}

template <typename CType>
void Current<CType>::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_id<<" {\n";
  msg_Indent();
#endif
  m_j.clear();
  Vertex_Vector::const_iterator vit(m_in.begin());
  // calculate outgoing momentum
  m_p=(*vit)->JA()->P()+(*vit)->JB()->P();
  // calculate subcurrents
  for (;vit!=m_in.end();++vit) 
    if (!(*vit)->Zero()) (*vit)->Evaluate();
  if (!m_out.empty() && !m_j.empty()) {
    AddPropagator();
    for (vit=m_out.begin();vit!=m_out.end();++vit)
      (*vit)->SetZero(m_j.empty());
  }
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
		 <<m_key<<"}("<<Flav()<<"){\n";
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
  return static_cast<Current<CType>*>(this)->AddJ(j); 
}

template <typename CType> const std::vector<CType> &Current_Base::J() const
{ 
  return static_cast<const Current<CType>*>(this)->J(); 
}

#include "C_Spinor.H"

template void Current_Base::AddJ(const CSpinor &j);
template const std::vector<CSpinor> &Current_Base::J() const;

template class Current<CSpinor>;

#include "C_Vector.H"

template void Current_Base::AddJ(const CVec4D &j);
template const std::vector<CVec4D> &Current_Base::J() const;

template class Current<CVec4D>;

#include "C_Tensor.H"

template void Current_Base::AddJ(const CAsT4D &j);
template const std::vector<CAsT4D> &Current_Base::J() const;

template class Current<CAsT4D>;
