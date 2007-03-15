#include "Vector_Current.H"

#include "Vertex.H"
#include "Tensor_Current.H"
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

static const double sqrttwo(sqrt(2.0));
static const double invsqrttwo(1.0/sqrttwo);

Vector_Current::Vector_Current(CDBG_Amplitude *const ampl,const Flavour &fl):
  Current(ampl,fl) {}

Vector_Current::~Vector_Current()
{
}

Vertex *Vector_Current::GetVertex(Current *const ja,Current *const jb)
{
  if ((ja->Type()==ct::vector && jb->Type()!=ct::spinor) ||
      (jb->Type()==ct::vector && ja->Type()!=ct::spinor)) {
    Vertex *v(new Vertex(this));
    v->SetJA(ja);
    v->SetJB(jb);
    return v;
  }
  return NULL;
}

void Vector_Current::SetCurrent()
{ 
#ifdef DEBUG__BG
  msg_Debugging()<<"set current "<<m_id<<" "<<p_ampl->JV()<<"\n";
#endif
  m_j.resize(1);
  m_j.back()=p_ampl->JV();
  for (Vertex_Vector::const_iterator vit=m_out.begin();
       vit!=m_out.end();++vit) (*vit)->SetZero(false);
}

void Vector_Current::AddJ(const CVec4D &j)
{
  PROFILE_HERE;
  for (CVec_Vector::iterator cit(m_j.begin());
       cit!=m_j.end();++cit)
    if (j(0)==(*cit)(0) && j(1)==(*cit)(1)) {
#ifdef DEBUG__BG
      msg_Debugging()<<"add "<<j<<"\n";
#endif
      *cit+=j;
      return;
    }
#ifdef DEBUG__BG
  msg_Debugging()<<"new "<<j<<"\n";
#endif
  m_j.push_back(j);
}

void Vector_Current::EvaluateTGC(Vector_Current *a,Vector_Current *b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<*a<<"(+)"<<*b<<"\n";
  msg_Indent();
#endif
  // triple gluon vertex
  for (CVec_Vector::const_iterator ait(a->m_j.begin());
       ait!=a->m_j.end();++ait) {
    for (CVec_Vector::const_iterator bit(b->m_j.begin());
	 bit!=b->m_j.end();++bit) {
      if ((*ait)(0)==(*bit)(1) || (*ait)(1)==(*bit)(0)) {
	// sum V_3(\pi_1,\pi_2) & V_3(\pi_2,\pi_1)
	CVec4D j(invsqrttwo*(*ait**bit)*(a->P()-b->P())
		 +sqrttwo*((*ait*b->P())**bit-(*bit*a->P())**ait));
	if ((*ait)(0)==(*bit)(1)) {
	  if ((*ait)(1)==(*bit)(0) && (*ait)(0)==(*ait)(1))
	    // color singlet, yields zero
	    continue;
	  // case 1: original ordering, V_3(\pi_1,\pi_2)
	  j(0)=(*bit)(0);
	  j(1)=(*ait)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"V case 1";
#endif
	  AddJ(j);
	}
	if ((*ait)(1)==(*bit)(0)) {
	  // case 2: reverse ordering, V_3(\pi_2,\pi_1)
	  j(0)=(*ait)(0);
	  j(1)=(*bit)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"V case 2";
#endif
	  AddJ(-j);
	}
      }
    }
  }
}

void Vector_Current::EvaluateQGC(Vector_Current *a,Tensor_Current *b,
				 const int mode)
{
#ifdef DEBUG__BG
  msg_Debugging()<<*a<<"(+)"<<*b<<" "<<mode<<"\n";
  msg_Indent();
#endif
  // gluon/tensor vertex
  for (CVec_Vector::const_iterator ait(a->m_j.begin());
       ait!=a->m_j.end();++ait) {
    for (CAsT_Vector::const_iterator bit(b->J().begin());
	 bit!=b->J().end();++bit) {
      if ((*ait)(0)==(*bit)(1) || (*ait)(1)==(*bit)(0)) {
	// sum V_4(\pi_1,\{\pi_2,\pi_3\}) & V_4(\{\pi_2,\pi_3\},\pi_1)
	CVec4D j(*ait**bit/2.0);
	if ((*ait)(0)==(*bit)(1)) {
	  if ((*ait)(1)==(*bit)(0) && (*ait)(0)==(*ait)(1))
	    // color singlet, yields zero
	    continue;
	  // case 1: original ordering, V_4(\pi_1,\{\pi_2,\pi_3\})
	  j(0)=(*bit)(0);
	  j(1)=(*ait)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"VT case 1";
#endif
	  AddJ(j);
	}
	if ((*ait)(1)==(*bit)(0)) {
	  // case 2: reverse ordering, V_4(\{\pi_2,\pi_3\},\pi_1)
	  j(0)=(*ait)(0);
	  j(1)=(*bit)(1);
#ifdef DEBUG__BG
	  msg_Debugging()<<"VT case 2";
#endif
	  AddJ(-j);
	}
      }
    }
  }
}

void Vector_Current::Evaluate(Vertex *const v)
{
  if (v->JA()->Type()==ct::vector) {
    if (v->JB()->Type()==ct::vector)
      EvaluateTGC((Vector_Current*)v->JA(),(Vector_Current*)v->JB());
    else if (v->JB()->Type()==ct::tensor)
      EvaluateQGC((Vector_Current*)v->JA(),(Tensor_Current*)v->JB(),1);
  }
  else if (v->JA()->Type()==ct::tensor) {
    if (v->JB()->Type()==ct::vector)
      EvaluateQGC((Vector_Current*)v->JB(),(Tensor_Current*)v->JA(),-1);
  }
}

void Vector_Current::Evaluate()
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
  if (!m_out.empty()) {
    // add propagator for off-shell leg
    double prop(1.0/m_p.Abs2());
#ifdef DEBUG__BG
    msg_Debugging()<<"propagator: "<<prop<<"\n";
#endif
    for (CVec_Vector::iterator jit(m_j.begin());
	 jit!=m_j.end();++jit) *jit*=prop;
    for (vit=m_out.begin();vit!=m_out.end();++vit)
      (*vit)->SetZero(m_j.empty());
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
  Print();
#endif
}

Complex Vector_Current::operator*(const Current &c) const
{
  PROFILE_HERE;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  if (c.Type()!=ct::vector) THROW(fatal_error,"Invalid current type.");
  Complex res(0.0);
  for (CVec_Vector::const_iterator jit1(m_j.begin());
       jit1!=m_j.end();++jit1)
    for (CVec_Vector::const_iterator jit2(((Vector_Current&)c).m_j.begin());
	 jit2!=((Vector_Current&)c).m_j.end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0)) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

void Vector_Current::Print() const
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
      msg_Debugging()<<m_j[i]<<" -> cc "<<m_j[i]*m_p<<"\n";
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

