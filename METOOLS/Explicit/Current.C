#include "METOOLS/Explicit/Current.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Current_Key
#define OBJECT_TYPE METOOLS::Current
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;
using namespace ATOOLS;

char METOOLS::ParticleType(const Flavour &fl)
{
  switch(fl.IntSpin()) {
  case 0: return 'S';
  case 1: return 'F';
  case 2: return 'V';
  case 4: return fl.IsDummy()?'P':'T';
  }
  msg_Error()<<METHOD<<"(): "<<fl<<std::endl;
  THROW(fatal_error,"Representation not implemented");
}

Current::Current(const Current_Key &key):
  m_fl(key.m_fl), m_key(0), m_order(2,0), m_cid(0), m_ntc(0),
  m_mass(m_fl.Mass()), m_width(m_fl.Width()), 
  m_msv(!IsZero(m_mass)), m_zero(true),
  m_dir(0), m_cut(0), m_osd(0), p_sub(NULL) {}

Current::~Current()
{
  ResetJ();
  for (Vertex_Vector::const_iterator vit(m_out.begin());
       vit!=m_out.end();++vit) (*vit)->ClearJ();
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) delete *vit;
}

void Current::FindPermutations()
{
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) (*vit)->FindPermutation();
}

void Current::InitPols(const Int_Vector &pols)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<pols<<") {\n";
#endif
  {
    msg_Indent();
    m_h.Init(pols);
    m_j.resize(m_h.N());
    for (Vertex_Vector::const_iterator vit(m_in.begin());
	 vit!=m_in.end();++vit) (*vit)->InitPols();
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

void Current::SetId(const Int_Vector &id)
{
  m_id=id;
  m_cid=0;
  for (size_t i(0);i<m_id.size();++i) m_cid+=1<<m_id[i];
}

std::string Current::PSInfo() const
{
  if (m_psinfo!="") return m_psinfo;
  std::string idt;
  for (size_t i(0);i<m_id.size();++i) m_psinfo+=ToString(m_id[i]);
  if (m_order[1]>0) m_psinfo+="_W"+ToString(m_order[1]);
  if (m_order[0]>0) m_psinfo+="_S"+ToString(m_order[0]);
  if (m_ntc>0) m_psinfo+="_T"+ToString(m_ntc);
  if (m_mass==0.0 && m_width==0.0) return m_psinfo;
  return m_psinfo+="["+ToString(m_mass)+","+ToString(m_width)+"]";
}

std::string Current::CLabel() const
{
  return "plain";
}

void Current::CollectGraphs
(Graph_Node *graph,const std::string &lastv) const
{
  if ((*graph)->empty()) {
    if (m_in.empty()) {
      graph->front()+=",j_"+ToString(CId());
      if (m_fl.IsAnti())
	graph->push_back("    \\fmf{"+CLabel()+"}{j_"
			 +ToString(CId())+","+lastv+"}");
      else 
	graph->push_back("    \\fmf{"+CLabel()+"}{"
			 +lastv+",j_"+ToString(CId())+"}");
    }
    else {
      for (Vertex_Vector::const_iterator vit(m_in.begin());
	   vit!=m_in.end();++vit) {
	Graph_Node *ngraph(new Graph_Node("",true));
	ngraph->pop_back();
	for (size_t j(0);j<graph->size();++j)
	  ngraph->push_back((*graph)[j]);
	if (m_fl.IsAnti()) 
	  ngraph->push_back("    \\fmf{"+CLabel()+"}{"
			    +(*vit)->VId()+","+lastv+"}");
	else
	  ngraph->push_back("    \\fmf{"+CLabel()+"}{"
			    +lastv+","+(*vit)->VId()+"}");
	(*vit)->CollectGraphs(ngraph);
	(*graph)().push_back(ngraph);
      }
    }
  }
  else {
    for (size_t i(0);i<(*graph)->size();++i)
      CollectGraphs((*graph)()[i],lastv);    
  }
}

void Current::CollectGraphs(Graph_Node *graph) const
{
  std::string lastv;
  for (String_Vector::reverse_iterator vit(graph->rbegin());
       vit!=graph->rend();++vit) {
    size_t pos(vit->rfind("%%"));
    if (pos!=std::string::npos) {
      lastv=vit->substr(pos+3);
      break;
    }
  }
  CollectGraphs(graph,lastv);
}

void Current::DetachOut(Vertex *const v)
{
  for (Vertex_Vector::iterator vit(m_out.begin());
       vit!=m_out.end();++vit)
    if (*vit==v) {
      m_out.erase(vit);
      return;
    }
  msg_Error()<<METHOD<<"(): Vertex '"<<v
	     <<"' not attached to current '"<<this<<"'"<<std::endl;
}

void Current::ResetJ()
{
  for (CObject_Matrix::iterator
	 ccit(m_j.begin());ccit!=m_j.end();++ccit) {
    for (CObject_Vector::const_iterator cit(ccit->begin());
	 cit!=ccit->end();++cit) (*cit)->Delete();
    ccit->clear();
  }
  m_zero=true;
}

void Current::Evaluate()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_id<<" {\n";
  msg_Indent();
#endif
  ResetJ();
  Vertex_Vector::const_iterator vit(m_in.begin());
  if (p_sub==NULL || m_id.size()>
      (p_sub->Sub()->In().front()->Info()->Mode()==1?2:1)) {
    // calculate outgoing momentum
    m_p=Vec4D();
    for (Current_Vector::const_iterator jit((*vit)->J().begin());
	 jit!=(*vit)->J().end();++jit) m_p+=(*jit)->P();
  }
  // calculate subcurrents
  for (;vit!=m_in.end();++vit) (*vit)->Evaluate();
  if (!m_out.empty() && !m_zero &&
      (p_sub==NULL || p_sub->Sub()!=this)) AddPropagator();
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
  Print();
#endif
}

void Current::ResetZero()
{
  for (Vertex_Vector::const_iterator vit(m_out.begin());
       vit!=m_out.end();++vit) if (!(*vit)->Zero()) return;
  ResetJ();
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) (*vit)->SetZero();
}

void Current::AddJ(CObject *const j)
{ 
  size_t i(j->H());
  for (CObject_Vector::iterator cit(m_j[i].begin());
       cit!=m_j[i].end();++cit)
    if (*j==**cit) { 
      (*cit)->Add(j);
      j->Delete();
      return; 
    }
  m_j[i].push_back(j);
  m_zero=false;
}

std::string Current_Key::Type() const
{
  return std::string(1,ParticleType(m_fl));
}

void Current::Print() const
{
  if (!msg_LevelIsDebugging()) return;
  std::string id(m_id.empty()?"<no entry>":ToString(m_id.front()));
  for (size_t i(1);i<m_id.size();++i) id+=","+ToString(m_id[i]);
  msg_Debugging()<<'['<<id<<"]"<<m_fid<<"{"<<m_id.size()<<","
		 <<m_key<<"}("<<m_order<<"|"<<m_ntc<<")("
		 <<(m_dir<=0?Flav():Flav().Bar())<<")"
		 <<(m_dir==0?"":m_dir>0?"I":"O")<<(m_cut?"c":"")
		 <<(p_sub?"S["+ToString(p_sub->Id())+
		    ToString(p_sub->Sub()->Id())+"]":"")<<"{\n";
  if (m_p!=Vec4D()) msg_Debugging()<<"m_p  : "<<m_p<<"\n";
  msg_Debugging()<<"m_j  :\n";
  for (size_t j(0);j<m_j.size();++j) {
    msg_Indent();
      for (size_t i(0);i<m_j[j].size();++i) 
	msg_Debugging()<<Format(m_j[j][i])<<"\n";
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

std::ostream &METOOLS::operator<<(std::ostream &str,const Current &c)
{
  return str<<'('<<c.Type()<<','<<c.Flav()<<','<<c.Id()
	    <<','<<c.FId()<<','<<c.Cut()<<')';
}

