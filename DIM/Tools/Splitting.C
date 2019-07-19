#include "DIM/Tools/Splitting.H"

#include "DIM/Tools/Amplitude.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>

using namespace DIM;
using namespace ATOOLS;

void Splitting::SetType()
{
  m_type=(p_c->Beam()?1:0)|(p_s->Beam()?2:0);
}

void Splitting::SetKinSpect(const Parton &p)
{
// p_s = spectator, p_c = emitter, p_n = emitted
  if   (p_c->Flav().IsbQuark() && p_s->Flav().IsbbarQuark())
    p_kinspec = (*p.Ampl())[2];
  else if(p_c->Flav().IsbbarQuark() && p_s->Flav().IsbQuark())
    p_kinspec = (*p.Ampl())[3];
  else p_kinspec = NULL;
}

void Splitting::SetKinVar()
{
  const double mw2            = sqr(Flavour(24).Mass());
  const ATOOLS::Vec4D pai     = p_c->Mom();
  const ATOOLS::Vec4D pb      = p_s->Mom();
  const ATOOLS::Vec4D pwtilde = p_kinspec->Mom();
  m_paipb   = pai*pb;
  m_Qprime2 = 2.*pai*pwtilde;
  const ATOOLS::Vec4D pminus  = pwtilde - mw2/m_Qprime2*pai;
  m_alpha   = pb*pminus / (pai*pminus);
}

bool Splitting::Allowed() const
{
  return std::find(p_s->S().begin(),p_s->S().end(),p_c)!=p_s->S().end();
}

namespace DIM {

  std::ostream &operator<<(std::ostream &s,const Splitting &p)
  {
    s<<"["<<(p.p_c?p.p_c->Id():0)<<"<->"<<(p.p_s?p.p_s->Id():0)
     <<"](kfac="<<p.m_kfac
     <<"){t="<<p.m_t<<",z="<<p.m_z<<",phi="<<p.m_phi<<"|";
    if (p.m_ci.size()) s<<p.m_ci[0]<<p.m_cj[0];
    for (size_t i(1);i<p.m_ci.size();++i) s<<";"<<p.m_ci[i]<<p.m_cj[i];
    return s<<"|"<<p.m_h[0]<<","<<p.m_h[1]<<","<<p.m_h[2]<<"}";
  }

}
