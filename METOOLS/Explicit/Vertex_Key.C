#include "METOOLS/Explicit/Vertex_Key.H"

#include "METOOLS/Explicit/Current.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace METOOLS;
using namespace ATOOLS;

AutoDelete_Vector<Vertex_Key> Vertex_Key::s_objects;

Vertex_Key::Vertex_Key
(const std::vector<Current*> &j,
 Current *const c,MODEL::Model_Base *const model,
 MODEL::Single_Vertex *const mv,const std::string &p,
 Vertex *const v,Color_Calculator *const cc,
 Lorentz_Calculator *const lc): 
  m_j(j), p_c(c), p_k(NULL), p_kt(NULL),
  p_model(model), p_mv(mv),
  m_p(p), m_n(0), m_d(0), p_v(v),
  p_cc(cc), p_lc(lc), p_dinfo(NULL)
{
}

Vertex_Key *Vertex_Key::New
(const std::vector<Current*> &j,
 Current *const c,MODEL::Model_Base *const model,
 MODEL::Single_Vertex *const mv,const std::string &p,
 Vertex *const v,Color_Calculator *const cc,
 Lorentz_Calculator *const lc)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new Vertex_Key(j,c,model,mv,p,v,cc,lc);
  Vertex_Key *k(s_objects.back());
  s_objects.pop_back();
  k->m_j=j;
  k->p_c=c;
  k->p_k=k->p_kt=NULL;
  k->p_model=model;
  k->p_mv=mv;
  k->m_p=p;
  k->m_n=0;
  k->m_d=0;
  k->p_v=v;
  k->p_cc=cc;
  k->p_lc=lc;
  k->p_dinfo=NULL;
  return k;
}

void Vertex_Key::Delete()
{
#ifdef USING__Threading
  s_objects.push_back(this);
#else
  delete this;
#endif
}

std::string Vertex_Key::Type() const
{
  std::string estr;
  for (size_t i(0);i<m_j.size();++i) estr+=m_j[i]->Type();
  return estr+p_c->Type();
}

const std::string &Vertex_Key::ID() const
{
  m_id.clear();
  for (size_t i(0);i<m_j.size();++i)
    m_id+=(m_j[i]?m_j[i]->Flav().IDName():
	   Flavour(p_dinfo->Type()?kf_photon:kf_gluon).IDName())+"|";
  if (p_c!=NULL) m_id+=p_c->Flav().Bar().IDName();
  return m_id;
}

ATOOLS::Flavour Vertex_Key::Fl(const size_t &i) const
{
  return m_j[i]?m_j[i]->Flav():Flavour(p_dinfo->Type()?kf_photon:kf_gluon);
}

Current *Vertex_Key::J(const size_t &i) const
{
  return m_j[i];
}
