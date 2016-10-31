#include "DIRE/Shower/Lorentz.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE DIRE::Kernel_Key
#define OBJECT_TYPE DIRE::Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "DIRE/Shower/Lorentz_FF.H"
#include "DIRE/Shower/Lorentz_II.H"

#include "DIRE/Tools/Parton.H"

using namespace DIRE;
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

void Lorentz::SetParams(Splitting &s,const PHASIC::Kin_Args &ff) const
{
  s.m_y=ff.m_y;
  s.m_x=ff.m_z;
  s.m_phi=ff.m_phi;
  s.m_mij2=p_ms->Mass2(m_fl[0]);
  s.m_mi2=p_ms->Mass2(m_fl[1]);
  s.m_mj2=p_ms->Mass2(m_fl[2]);
  s.m_mk2=p_ms->Mass2(s.p_s->Flav());
  s.m_Q2=dabs((s.p_c->Mom()+s.p_n->Mom()+s.p_s->Mom()).Abs2()
	      -s.m_mi2-s.m_mj2-s.m_mk2);
  s.p_sk=p_sk;
}

int Lorentz::Update(Splitting &s,const PHASIC::Kin_Args &ff,
		    const int mode,const Vec4D &pl) const
{
  if (ff.m_lam.size())
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      (*s.p_c->Ampl())[i]->SetMom
	(ff.m_lam*(*s.p_c->Ampl())[i]->Mom());
  ATOOLS::Vec4D pc(s.p_c->Mom()), ps(s.p_s->Mom());
  if (s.p_c->Out(0)==NULL) s.p_c->SetFlav(m_fl[1]);
  s.p_c->SetMom(ff.m_pi);
  s.p_s->SetMom(ff.m_pk);
  if (s.p_n==NULL) {
    s.p_n = new Parton(s.p_c->Ampl(),m_fl[2],ff.m_pj);
    s.p_n->SetId(s.p_n->Counter());
    s.p_c->Ampl()->Add(s.p_n);
    if (m_fl.size()>3) {
      s.p_l = new Parton(s.p_c->Ampl(),m_fl[3],pl);
      s.p_l->SetId(s.p_l->Counter());
      s.p_c->Ampl()->Add(s.p_l);
    }
  }
  else {
    if (s.p_n->Out(0)==NULL) s.p_n->SetFlav(m_fl[2]);
    s.p_n->SetMom(ff.m_pj);
  }
  if (mode&1) return 1;
  int stat(s.p_c->Out(0)==NULL);
  if (stat) {
    Jet_Finder *jf=s.p_c->Ampl()->JF<Jet_Finder>();
    if (jf) {
      Cluster_Amplitude *ampl(s.p_c->Ampl()->GetAmplitude());
      if (jf->JC()->Jets(ampl)) stat=0;
      if (stat) s.p_c->Ampl()->SetJF(NULL);
      msg_Debugging()<<(stat?"no ":"")<<"jet veto\n";
      ampl->Delete();
    }
  }
  return stat;
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
  s.m_Q2=dabs((s.p_c->Mom()+s.p_s->Mom()).Abs2()
	      -s.m_mi2-s.m_ml2-s.m_mj2-s.m_mk2);
  s.m_eta=s.p_c->GetXB();
  return true;
}
