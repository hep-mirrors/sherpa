#include "METOOLS/Explicit/Vertex_Key.H"

#include "METOOLS/Explicit/Current.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;

Vertex_Key::Vertex_Key
(Current *const a,Current *const b,Current *const e,
 Current *const c,MODEL::Model_Base *const model,
 MODEL::Single_Vertex *const mv,const std::string &p,
 Vertex *const v,Color_Calculator *const cc,
 Lorentz_Calculator *const lc): 
  p_a(a), p_b(b), p_c(c), p_e(e), p_k(NULL), p_kt(NULL),
  p_model(model), p_mv(mv), m_p(p), m_n(0), p_v(v),
  p_cc(cc), p_lc(lc), p_dinfo(NULL)
{
}

std::string Vertex_Key::Type() const
{
  std::string estr;
  if (p_e!=NULL) estr=p_e->Type();
  return std::string(1,p_a->Type())+p_b->Type()+estr+p_c->Type();
}

ATOOLS::Flavour Vertex_Key::FlA() const
{
  if (p_a) return p_a->Flav();
  if (p_dinfo) return Flavour(kf_gluon);
  return Flavour(kf_none);
}

ATOOLS::Flavour Vertex_Key::FlB() const
{
  if (p_b) return p_b->Flav();
  if (p_dinfo) return Flavour(kf_gluon);
  return Flavour(kf_none);
}

std::string Vertex_Key::ID(const int dir) const
{
  Flavour a(FlA()), b(FlB());
  std::string estr, cstr;
  if (p_e!=NULL) estr="{"+p_e->Flav().IDName()+"}";
  if (p_c!=NULL) cstr="{"+p_c->Flav().IDName()+"}";
  return '{'+a.IDName()+"}{"+b.IDName()+'}'+estr+cstr;
}

Vertex_Key Vertex_Key::SwapAB() const
{ 
  return Vertex_Key(p_b,p_a,p_e,p_c,p_model,p_mv,m_p,p_v,p_cc,p_lc); 
}

Vertex_Key Vertex_Key::SwapBE() const
{ 
  return Vertex_Key(p_a,p_e,p_b,p_c,p_model,p_mv,m_p,p_v,p_cc,p_lc); 
}

Vertex_Key Vertex_Key::SwapEA() const
{ 
  return Vertex_Key(p_e,p_b,p_a,p_c,p_model,p_mv,m_p,p_v,p_cc,p_lc); 
}

Current *Vertex_Key::J(const size_t &i) const
{
  switch (i) {
  case 0: return p_a;
  case 1: return p_b;
  case 2: return p_e;
  default: THROW(fatal_error,"Invalid index "+ToString(i));
  }
  return NULL;
}

bool Vertex_Key::operator<(const Vertex_Key &k) const
{
  if (!p_e&&k.p_e) return true;
  if (p_e&&!k.p_e) return false;
  if (p_a<k.p_a) return true;
  if (p_a>k.p_a) return false;
  if (p_b<k.p_b) return true;
  if (p_b>k.p_b) return false;
  if (p_e) {
    if (p_e<k.p_e) return true;
    if (p_e>k.p_e) return false;
  }
  return false;
}

