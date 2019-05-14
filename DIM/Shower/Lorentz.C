#include "DIM/Shower/Lorentz.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE DIM::Kernel_Key
#define OBJECT_TYPE DIM::Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Single_Vertex.H"
#include "DIM/Shower/Shower.H"
#include "DIM/Shower/Lorentz_FF.H"
#include "DIM/Shower/Lorentz_II.H"
#include "DIM/Tools/Parton.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace DIM;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Lorentz::Lorentz(const Kernel_Key &k,const int type):
  p_sk(k.p_k), m_type(type)
{
  m_fl[0]=k.p_v->in[0].Bar();
  if (k.m_mode==0) {
    m_fl[1]=k.p_v->in[1];
    m_fl[2]=k.p_v->in[2];
  }
  else {
    m_fl[1]=k.p_v->in[2];
    m_fl[2]=k.p_v->in[1];
  }

  Settings& s = Settings::GetMainSettings();
  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")    m_dipole_case = EXTAMP::DipoleCase::CS;
  else if (dipole_string == "IDa")   m_dipole_case = EXTAMP::DipoleCase::IDa;  // ee > bbWW with mapping a
  else if (dipole_string == "IDb")   m_dipole_case = EXTAMP::DipoleCase::IDb;  // ee > bbWW with mapping b
  else if (dipole_string == "IDin")  m_dipole_case = EXTAMP::DipoleCase::IDin; // pp > bbWW
  else if (dipole_string == "RES")   m_dipole_case = EXTAMP::DipoleCase::RES;  // ee > guu
  else if (dipole_string == "ID")    m_dipole_case = EXTAMP::DipoleCase::ID;   // ee > guu
  else                               m_dipole_case = EXTAMP::DipoleCase::CS;
}

Lorentz::~Lorentz()
{
}

double Lorentz::AsymmetryFactor(const Splitting &s) const
{
  return 1.0;
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
		    const int mode) const
{
  if (ff.m_lam.size())
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      (*s.p_c->Ampl())[i]->SetMom
	(ff.m_lam*(*s.p_c->Ampl())[i]->Mom());
  ATOOLS::Vec4D pc(s.p_c->Mom()), ps(s.p_s->Mom());
  s.p_c->SetFlav(m_fl[1]);
  s.p_c->SetMom(ff.m_pi);
  s.p_s->SetMom(ff.m_pk);
  Parton *n(NULL);
  if (s.p_n==NULL) {
    s.p_n = n = new Parton(s.p_c->Ampl(),m_fl[2],ff.m_pj);
    s.p_n->SetId(s.p_n->Counter());
    s.p_c->Ampl()->Add(n);
  }
  else {
    s.p_n->SetFlav(m_fl[2]);
    s.p_n->SetMom(ff.m_pj);
  }
  if (mode&1) return 1;
  int stat(1);
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
  s.m_mk2=p_ms->Mass2(s.p_s->Flav());
  s.m_Q2=dabs((s.p_c->Mom()+s.p_s->Mom()).Abs2()
	      -s.m_mi2-s.m_mj2-s.m_mk2);
  s.m_eta=s.p_c->GetXB();
  return true;
}
