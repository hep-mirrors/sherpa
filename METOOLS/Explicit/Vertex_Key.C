#include "METOOLS/Explicit/Vertex_Key.H"

#include "METOOLS/Explicit/Current.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;

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

std::string Vertex_Key::Type() const
{
  std::string estr;
  for (size_t i(0);i<m_j.size();++i) estr+=m_j[i]->Type();
  return estr+p_c->Type();
}

std::string Vertex_Key::ID() const
{
  std::string estr;
  for (size_t i(0);i<m_j.size();++i)
    estr+="{"+(m_j[i]?m_j[i]->Flav().IDName():
	       Flavour(p_dinfo->Type()?kf_photon:kf_gluon).IDName())+"}";
  if (p_c!=NULL) estr+="{"+p_c->Flav().Bar().IDName()+"}";
  return estr;
}

ATOOLS::Flavour Vertex_Key::Fl(const size_t &i) const
{
  return m_j[i]?m_j[i]->Flav():Flavour(p_dinfo->Type()?kf_photon:kf_gluon);
}

Current *Vertex_Key::J(const size_t &i) const
{
  return m_j[i];
}
