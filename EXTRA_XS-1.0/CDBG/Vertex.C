#include "Vertex.H"

#include "Message.H"
#include "Exception.H"
#include "STL_Tools.H"
#include "MyStrStream.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE EXTRAXS::Vertex_Key
#define OBJECT_TYPE EXTRAXS::Vertex
#include "Getter_Function.C"

using namespace EXTRAXS;
using namespace ATOOLS;

Vertex::Vertex(Current_Base *const c): 
  p_a(NULL), p_b(NULL), p_c(c), 
  m_zero(true), m_sign(false), m_fperm(0)
{
  if (c!=NULL) p_c->AttachIn(this);
}

Vertex::~Vertex()
{
}

void Vertex::FindPermutation()
{
  m_fperm=0;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Int_Vector id(p_c->Id()), fid(p_c->FId());
  Int_Vector pid(p_a->Id()), pfid(p_a->FId());
  if (1) {//p_b->Id().front()>id.front()) {
    pid.insert(pid.end(),p_b->Id().begin(),p_b->Id().end());
    pfid.insert(pfid.end(),p_b->FId().begin(),p_b->FId().end());
  }
  else {
    pid.insert(pid.begin(),p_b->Id().begin(),p_b->Id().end());
    pfid.insert(pfid.begin(),p_b->FId().begin(),p_b->FId().end());
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"  pid = "<<pid<<", pfid = "<<pfid<<"\n";
  msg_Debugging()<<"  id  = "<<id<<", fid  = "<<fid<<"\n";
#endif
  for (size_t i(0);i<id.size();++i) {
    for (size_t j(0);j<pid.size();++j) 
      if (pid[j]==id[i] && i!=j) {
	for (size_t k(j);k!=i;) {
	  size_t l(j>i?k-1:k+1);
 	  m_fperm+=pfid[k]==1&&pfid[l]==1;
#ifdef DEBUG__BG
	  if (pfid[k]==1 && pfid[l]==1)
	    msg_Debugging()<<"  swap "<<pid[l]<<" & "<<pid[k]<<"\n";
#endif
	  std::swap<int>(pid[k],pid[l]);
	  std::swap<int>(pfid[k],pfid[l]);
	  k=l;
	}
	break;
      }
  }
  m_sign=m_fperm%2==1;
#ifdef DEBUG__BG
  msg_Debugging()<<"} => "<<*this<<"\n";
#endif
}

std::string Vertex::VId() const
{
  return "v_"+ToString(p_a->CId())+"_"
    +ToString(p_b->CId())+"_"+ToString(p_c->CId());
}

std::string Vertex::VLabel() const
{
  return "decor.size=0ex,label=$\\white F="+ToString(m_fperm)+"$";
}

void Vertex::CollectGraphs(Graph_Node *graph) const
{
  graph->push_back("    \\fmfv{"+VLabel()+"}{"+VId()+"}");
  graph->push_back("    %% "+VId());
  p_a->CollectGraphs(graph);
  p_b->CollectGraphs(graph);
}

std::string Vertex_Key::Type() const
{
  return m_model+"_"+p_a->Type()+p_b->Type()+p_c->Type();
}

bool Vertex_Key::operator<(const Vertex_Key &k) const
{
  if (p_a<k.p_a) return true;
  if (p_a>k.p_a) return false;
  if (p_b<k.p_b) return true;
  if (p_b>k.p_b) return false;
  if (p_c<k.p_c) return true;
  if (p_c>k.p_c) return false;
  if (m_model<k.m_model) return true;
  return false;
}

std::ostream &EXTRAXS::operator<<(std::ostream &str,const Vertex &v)
{
  if (v.JA()!=NULL) 
    str<<'{'<<v.JA()->Type()<<','<<v.JA()->Flav()<<'}'<<v.JA()->Id();
  if (v.JB()!=NULL) 
    str<<"(+){"<<v.JB()->Type()<<','<<v.JB()->Flav()<<'}'<<v.JB()->Id();
  return str<<" {"<<v.FPerm()<<","<<v.Sign()<<"}";
}
