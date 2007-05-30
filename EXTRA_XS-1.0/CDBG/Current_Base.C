#include "Current.H"

#include "Vertex.H"
#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "STL_Tools.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE EXTRAXS::Current_Key
#define OBJECT_TYPE EXTRAXS::Current_Base
#include "Getter_Function.C"

using namespace EXTRAXS;
using namespace ATOOLS;

size_t EXTRAXS::MakeId(const Int_Vector &ids)
{
  size_t id(0);
  for (size_t i(0);i<ids.size();++i) 
    if (ids[i]>0) id+=1<<i;
  return id;
}

Int_Vector EXTRAXS::MakeId(const size_t &id,const size_t &n)
{
  size_t ic(id);
  Int_Vector ids(n,-1);
  for (size_t i(0);i<ids.size();++i) {
    size_t c(1<<i);
    if (ic&c) {
      ids[i]=1;
      ic-=c;
    }
  }
  if (ic!=0) THROW(fatal_error,"Invalid particle number");
  return ids;
}

char EXTRAXS::ParticleType(const Flavour &fl)
{
  switch(fl.IntSpin()) {
  case 1: return 'S';
  case 2: return 'V';
  case 4: return 'T';
  }
  THROW(fatal_error,"Representation not implemented");
}

Current_Base::Current_Base(const Flavour &fl):
  m_fl(fl), m_key(0) {}

Current_Base::~Current_Base()
{
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

std::string Current_Key::Type() const
{
  return m_model+"_"+ParticleType(m_fl);
}

std::ostream &EXTRAXS::operator<<(std::ostream &str,const Current_Base &c)
{
  return str<<'('<<c.Type()<<','<<c.Flav()<<','<<c.Id()<<','<<c.FId()<<')';
}

