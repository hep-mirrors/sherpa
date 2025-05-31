#include "MALARIC/Shower/Lorentz.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE MALARIC::Kernel_Key
#define OBJECT_TYPE MALARIC::Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Single_Vertex.H"
#include "MALARIC/Shower/Shower.H"
#include "MALARIC/Shower/Lorentz_FS.H"
#include "MALARIC/Shower/Lorentz_IS.H"
#include "MALARIC/Tools/Parton.H"

using namespace MALARIC;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Lorentz::Lorentz(const Kernel_Key &k,const int type):
  p_sk(k.p_k), m_type(type), m_fl(3)
{
  if (k.p_v==NULL) {
    m_fl=k.m_fl;
    m_fl[0]=m_fl[0].Bar();
    return;
  }
  m_fl[0]=k.p_v->in[0].Bar();
  if (k.m_mode==0) {
    m_fl[1]=k.p_v->in[1];
    m_fl[2]=k.p_v->in[2];
  }
  else {
    m_fl[1]=k.p_v->in[2];
    m_fl[2]=k.p_v->in[1];
  }
}

Lorentz::~Lorentz()
{
}

void Lorentz::SetMS(const ATOOLS::Mass_Selector *ms)
{
  p_ms=ms;
}

double Lorentz::AsymmetryFactor(const Splitting &s) const
{
  return 1.0;
}

void Lorentz::SetParams(Splitting &s) const
{
  s.m_mij2=p_ms->Mass2(m_fl[0]);
  s.m_mi2=p_ms->Mass2(m_fl[1]);
  s.m_mj2=p_ms->Mass2(m_fl[2]);
  s.m_mk2=p_ms->Mass2(s.p_s->Flav());
  s.m_q2=(s.p_c->Mom()+s.p_n->Mom()+s.p_s->Mom()).Abs2();
  s.m_Q2=dabs(s.m_q2-s.m_mi2-s.m_mj2-s.m_mk2);
  s.p_sk=p_sk;
}

int Lorentz::Update(Splitting &s,const int mode) const
{
  if (s.m_lam.size())
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      (*s.p_c->Ampl())[i]->SetMom
	(s.m_lam*(*s.p_c->Ampl())[i]->Mom());
  if (s.m_p.size())
    for (size_t i(0);i<s.m_p.size();++i)
      if (s.m_p[i]!=Vec4D())
	(*s.p_c->Ampl())[i]->SetMom(s.m_p[i]);
  ATOOLS::Vec4D pc(s.p_c->Mom()), ps(s.p_s->Mom());
  s.p_c->SetFlav(m_fl[1]);
  s.p_c->SetMom(s.m_pi);
  if (s.p_n==NULL) {
    s.p_n = new Parton(s.p_c->Ampl(),m_fl[2],s.m_pj);
    s.p_n->SetId(s.p_n->Counter());
    s.p_c->Ampl()->Add(s.p_n);
    if (m_fl.size()>3) {
      s.p_l = new Parton(s.p_c->Ampl(),m_fl[3],s.m_pl);
      s.p_l->SetId(s.p_l->Counter());
      s.p_c->Ampl()->Add(s.p_l);
    }
  }
  else {
    s.p_n->SetFlav(m_fl[2]);
    s.p_n->SetMom(s.m_pj);
  }
  if (mode&2) return 1;
  return 1;
}

bool Lorentz::Allowed(const Splitting &s) const
{
  return s.m_type==m_type && s.p_c->Flav()==m_fl[0];
}

bool Lorentz::SetLimits(Splitting &s) const
{
  s.m_t0=p_sk->PS()->TMin(s.m_type&1);
  s.m_mij2=p_ms->Mass2(m_fl[0]);
  s.m_mi2=p_ms->Mass2(m_fl[1]);
  s.m_mj2=p_ms->Mass2(m_fl[2]);
  if (m_fl.size()>3) s.m_ml2=p_ms->Mass2(m_fl[3]);
  s.m_mk2=p_ms->Mass2(s.p_s->Flav());
  s.m_q2=(s.p_c->Mom()+s.p_s->Mom()).Abs2();
  s.m_Q2=dabs(s.m_q2-s.m_mi2-s.m_ml2-s.m_mj2-s.m_mk2);
  if (!s.m_clu) {
    s.m_Kt=Vec4D();
    s.m_iink=0;
    const Amplitude &a(*s.p_c->Ampl());
    for (size_t i(0);i<a.size();++i)
      if (s.m_rcl[i]&2) {
	s.m_Kt+=a[i]->Mom();
	if (a[i]==s.p_c) s.m_iink=1;
      }
    if (s.m_iink) s.m_Kt=-s.m_Kt;
  }
  s.m_eta=s.p_c->GetXB();
  return true;
}
