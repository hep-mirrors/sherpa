#include "COMIX/Amplitude/Current_Base.H"

#include "COMIX/Amplitude/Vertex_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE COMIX::Current_Key
#define OBJECT_TYPE COMIX::Current_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace COMIX;
using namespace ATOOLS;

char COMIX::ParticleType(const Flavour &fl)
{
  switch(fl.IntSpin()) {
  case 0: return 'S';
  case 1: return 'F';
  case 2: return 'V';
  case 4: return 'T';
  }
  msg_Error()<<METHOD<<"(): "<<fl<<std::endl;
  THROW(fatal_error,"Representation not implemented");
}

Current_Base::Current_Base(const Current_Key &key):
  m_fl(key.m_fl), m_key(0), m_oew(0), m_oqcd(0), m_ntc(0),
  m_mass(m_fl.Mass()), m_width(m_fl.Width()), 
  m_msv(!IsZero(m_mass)), m_zero(true), m_dir(0),
  m_cut(0), m_osd(0) {}

Current_Base::~Current_Base()
{
  for (Vertex_Vector::const_iterator vit(m_out.begin());
       vit!=m_out.end();++vit)
    if ((*vit)->JA()==this) (*vit)->SetJA(NULL);
    else if ((*vit)->JB()==this) (*vit)->SetJB(NULL);
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) delete *vit;
}

void Current_Base::FindPermutations()
{
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) (*vit)->FindPermutation();
}

size_t Current_Base::CId() const
{
  size_t id(0);
  for (size_t i(0);i<m_id.size();++i) id+=1<<m_id[i];
  return id;
}

std::string Current_Base::PSInfo() const
{
  if (m_psinfo!="") return m_psinfo;
  std::string idt;
  for (size_t i(0);i<m_id.size();++i) m_psinfo+=ToString(m_id[i]);
  if (m_mass==0.0 && m_width==0.0) return m_psinfo;
  return m_psinfo+="["+ToString(m_mass)+","+ToString(m_width)+"]";
}

std::string Current_Base::CLabel() const
{
  return "plain";
}

void Current_Base::CollectGraphs
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

void Current_Base::CollectGraphs(Graph_Node *graph) const
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

void Current_Base::DetachOut(Vertex_Base *const v)
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

std::string Current_Key::Type() const
{
  return std::string(1,ParticleType(m_fl));
}

std::ostream &COMIX::operator<<(std::ostream &str,const Current_Base &c)
{
  return str<<'('<<c.Type()<<','<<c.Flav()<<','<<c.Id()
	    <<','<<c.FId()<<','<<c.Cut()<<')';
}

