#include "Tensor_Current.H"

#include "Vertex.H"
#include "CDBG_Amplitude.H"
#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

// #define DEBUG__BG

using namespace EXTRAXS;
using namespace ATOOLS;

Tensor_Current::Tensor_Current(CDBG_Amplitude *const ampl,const Flavour &fl):
  Current(ampl,fl) {}

Tensor_Current::~Tensor_Current()
{
}

Vertex *Tensor_Current::GetVertex(Current *const ja,Current *const jb)
{
  if (ja->Type()==ct::vector && jb->Type()==ct::vector) {
    Vertex *v(new Vertex(this));
    v->SetJA(ja);
    v->SetJB(jb);
    return v;
  }
  return NULL;
}

void Tensor_Current::SetCurrent()
{ 
  THROW(fatal_error,"External tensor particle.");
}

void Tensor_Current::AddJ(const CAsT4D &j)
{
  PROFILE_HERE;
  bool found(false);
  for (CAsT_Vector::iterator cit(m_j.begin());
       cit!=m_j.end();++cit)
    if (j(0)==(*cit)(0) && j(1)==(*cit)(1)) {
      found=true;
#ifdef DEBUG__BG
      msg_Debugging()<<"add "<<j<<"\n";
#endif
      *cit+=j;
    }
  if (!found) {
#ifdef DEBUG__BG
    msg_Debugging()<<"new "<<j<<"\n";
#endif
    m_j.push_back(j);
  }
}

void Tensor_Current::EvaluateQGC(Vector_Current *a,Vector_Current *b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<a->Id()<<"(+)"<<b->Id()<<"\n";
  msg_Indent();
#endif
  for (CVec_Vector::const_iterator jit1(a->J().begin());
       jit1!=a->J().end();++jit1) 
    for (CVec_Vector::const_iterator jit2(b->J().begin());
	 jit2!=b->J().end();++jit2) {
      // sum T_3(\pi_1,\pi_2) & T_3(\pi_2,\pi_1)
      if ((*jit1)(0)==(*jit2)(1)) {
	if ((*jit1)(1)==(*jit2)(0) && (*jit1)(0)==(*jit1)(1))
	  // color singlet, yields zero
	  continue;
#ifdef DEBUG__BG
	msg_Debugging()<<"TC case 1";
#endif
	AddJ(CAsT4D(*jit1,*jit2));
      }
      if ((*jit1)(1)==(*jit2)(0)) {
#ifdef DEBUG__BG
	msg_Debugging()<<"TC case 2";
#endif
	AddJ(CAsT4D(*jit2,*jit1));
      }
    }
}

void Tensor_Current::Evaluate(Vertex *const v)
{
  EvaluateQGC((Vector_Current*)v->JA(),(Vector_Current*)v->JB());
}

void Tensor_Current::Evaluate()
{
  PROFILE_HERE;
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
    if (!(*vit)->Zero()) Evaluate(*vit);
  for (vit=m_out.begin();vit!=m_out.end();++vit)
    (*vit)->SetZero(m_j.empty());
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
  Print();
#endif
}

Complex Tensor_Current::operator*(const Current &c) const
{
  THROW(fatal_error,"Multiplication of tensor particles not allowed");
  return 0.0;
}

void Tensor_Current::Print() const
{
  if (!msg.LevelIsDebugging()) return;
  std::string id(m_id.empty()?"<no entry>":ToString(m_id.front()));
  for (size_t i(1);i<m_id.size();++i) id+=","+ToString(m_id[i]);
  msg_Debugging()<<'['<<id<<"]("<<m_id.size()<<","
		 <<m_key<<")("<<Flav()<<"){\n";
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

