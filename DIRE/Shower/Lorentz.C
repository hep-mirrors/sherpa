#include "DIRE/Shower/Lorentz.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE DIRE::Kernel_Key
#define OBJECT_TYPE DIRE::Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "DIRE/Shower/Lorentz_FF.H"
#include "DIRE/Shower/Lorentz_II.H"
#include "DIRE/Tools/Parton.H"
#include "ATOOLS/Org/Scoped_Settings.H"

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

  Settings& s = Settings::GetMainSettings();
  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")    m_dipole_case = EXTAMP::DipoleCase::CS;
  else if (dipole_string == "IDa")   m_dipole_case = EXTAMP::DipoleCase::IDa;  // ee > bbWW with mapping a
  else if (dipole_string == "IDb")   m_dipole_case = EXTAMP::DipoleCase::IDb;  // ee > bbWW with mapping b
  else if (dipole_string == "IDin")  m_dipole_case = EXTAMP::DipoleCase::IDin; // pp > bbWW
  else if (dipole_string == "RES")   m_dipole_case = EXTAMP::DipoleCase::RES;  // ee > guu
  else if (dipole_string == "ID")    m_dipole_case = EXTAMP::DipoleCase::ID;   // ee > guu
  else                               m_dipole_case = EXTAMP::DipoleCase::CS;

  m_t_cutoff = s["CSS_FS_PT2MIN"].Get<double>();
  m_evol     = s["CSS_EVOLUTION_SCHEME"].Get<int>();
  m_maxem    = s["CSS_MAXEM"].Get<size_t>();
}

Lorentz::~Lorentz()
{
}

void Lorentz::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
}

void Lorentz::SetParams(Splitting &s,const PHASIC::Kin_Args &ff) const
{
  s.m_ff=ff;
  s.m_y=ff.m_y;
  s.m_x=ff.m_z;
  s.m_phi=ff.m_phi;
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
  ATOOLS::Vec4D pc(s.p_c->Mom()), ps(s.p_s->Mom());
  if (s.p_c->Out(0)==NULL) s.p_c->SetFlav(m_fl[1]);
  switch(m_dipole_case){
    case EXTAMP::IDa:
    {
      s.p_c->SetMom(s.m_pk);
      s.p_kinspec->SetMom(s.m_pi);
      break;
    }
    default:
    {
      s.p_c->SetMom(s.m_pi);
      s.p_s->SetMom(s.m_pk);
    }
  }
  if (s.p_n==NULL) {
    s.p_n = new Parton(s.p_c->Ampl(),m_fl[2],s.m_pj);
    if(m_dipole_case == EXTAMP::IDa) s.p_n->SetKinSpectID(s.p_c->GetKinSpectID());
    s.p_n->SetId(s.p_n->Counter());
    s.p_c->Ampl()->Add(s.p_n);
    if (m_fl.size()>3) {
      s.p_l = new Parton(s.p_c->Ampl(),m_fl[3],s.m_pl);
      s.p_l->SetId(s.p_l->Counter());
      s.p_c->Ampl()->Add(s.p_l);
    }
  }
  else {
    if (s.p_n->Out(0)==NULL) s.p_n->SetFlav(m_fl[2]);
    s.p_n->SetMom(s.m_pj);
  }
  if (mode&2) return 1;
  int stat(s.p_c->Out(0)==NULL);
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
  s.m_q2=(s.p_c->Mom()+s.p_s->Mom()).Abs2();
  s.m_Q2=dabs(s.m_q2-s.m_mi2-s.m_ml2-s.m_mj2-s.m_mk2);
  s.m_eta=s.p_c->GetXB();

  if(m_dipole_case==EXTAMP::IDa){
    /* if shower is invoked, but not executed -> happens in matching, but does to cause trouble */
    if(m_maxem==0) { return true; }

    if(m_evol==1){
      const double k02     = m_t_cutoff;
      const double mw2     = sqr(Flavour(24).Mass());
      const double paipb   = s.m_paipb;
      const double Qprime2 = s.m_Qprime2;
      const double alpha   = s.m_alpha;
      s.m_zmax   = (-2.*k02*paipb*Qprime2 + paipb*sqr(Qprime2) - alpha*pow(Qprime2,3.) +
                   sqrt(alpha*pow(Qprime2,4.)*(4.*k02*paipb + alpha*sqr(Qprime2)))) /
                   (paipb*sqr(Qprime2));
    }
    else if(m_evol==2){
      const double tmin     = m_t_cutoff;
      const double mw2     = sqr(Flavour(24).Mass());
      const double Qprime2 = s.m_Qprime2;
      const double radicand = sqr(1.+tmin/Qprime2)/4. - tmin/Qprime2*(1+mw2/Qprime2);

      if(radicand >= 0) s.m_zmax = (1.+tmin/Qprime2)/2 + sqrt(radicand);
      else              s.m_zmax = (1.+tmin/Qprime2)/2;

      /* following line is rarely relevant: only if Qprime2 is extremely small */
      if(s.m_zmax>1.) s.m_zmax = 0.99999;        // TODO
    }
    /* in case zmax is slightly above 1, due to numerics */
    if(s.m_zmax>1. && IsEqual(s.m_zmax,1.,1.e-8))           s.m_zmax=0.99999;
    if(!(s.m_zmax>0.) || !(s.m_zmax<=1.))    THROW(fatal_error, "zmax wrong");
  }
  return true;
}
