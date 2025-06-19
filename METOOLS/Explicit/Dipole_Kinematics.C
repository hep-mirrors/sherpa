#include "METOOLS/Explicit/Dipole_Kinematics.H"

#include "METOOLS/Explicit/Vertex.H"
#include "PDF/Main/NLOMC_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PHASIC++/Channels/Antenna_Kinematics.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"

using namespace METOOLS;
using namespace ATOOLS;

Dipole_Kinematics::Dipole_Kinematics
(Dipole_Info *const info,Current *const i,Current *const j,
 Current *const k,Current *const ijt,Current *const kt):
  p_i(i), p_j(j), p_k(k), p_ijt(ijt), p_kt(kt),
  m_type(0), m_swap(0), m_trig(1), m_z(0.), m_y(0.), m_kt2(0.), m_Q2(0.),
  p_info(info), m_mi2(0.0), m_mj2(0.0),
  m_mij2(sqr(p_ijt->Flav().Mass())), m_mk2(sqr(p_k->Flav().Mass())),
  m_ym(0.0), m_yp(1.0), m_f(0.0), m_a(0.0),
  p_subevt(nullptr), p_nlomc(nullptr),
  p_softrecoil(nullptr), p_collrecoil(nullptr)
{
  if (p_i) m_mi2=sqr(p_i->Flav().Mass());
  if (p_j) m_mj2=sqr(p_j->Flav().Mass());
  m_phase[1]=m_phase[0]=0.0;
  m_res[2]=m_res[1]=m_res[0]=0.0;
  if ((p_i && p_i->Direction()==0) ||
      (p_j && p_j->Direction()==0) ||
      p_k->Direction()==0)
    THROW(fatal_error,"Missing current information");
  if (p_k->Direction()>0) m_type|=2;
  if ((p_i && p_i->Direction()>0) ||
      (p_j && p_j->Direction()>0)) {
    if (p_i && p_j && p_j->Direction()>0) {
      std::swap<Current*>(p_i,p_j);
      m_swap=1;
    }
    m_type|=1;
  }
}

Dipole_Kinematics::~Dipole_Kinematics()
{
}

void Dipole_Kinematics::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  p_info->SetSubType(p_nlomc->SubtractionType());
  if (p_nlomc->SubtractionType()==1) p_info->SetKappa(1.0);
}

double Dipole_Kinematics::Lam
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

void Dipole_Kinematics::Evaluate()
{
  m_pi=p_i->P();
  m_pj=p_j->P();
  m_pk=p_k->P();
  m_Q=m_pi+m_pj+m_pk;
  m_Q2=m_Q.Abs2();
  if (p_info->SubType()==subscheme::code::Alaric)
    EvaluateAlaricKinematics();
  else
    EvaluateCSSKinematics();
  m_trig = m_y/m_yp < ((p_nlomc) ? 1.0 : p_info->AMax(m_type));
  if (m_trig) m_trig=m_kt2<p_info->KT2Max();
  if (p_info->AMin()>0.0) {
    if (m_y<p_info->AMin()) p_info->SetStat(0);
  }
  else {
    if (m_kt2<-p_info->AMin()) p_info->SetStat(0);
  }
  if (p_info->Stat() &&
      !IsEqual(p_i->P()+p_j->P()+p_k->P(),p_ijt->P()+p_kt->P())) {
    msg_Error()<<METHOD<<"(): Momentum not conserved in type "<<m_cur.back()->SubType()
	       <<" {\n  before "<<p_i->P()+p_j->P()+p_k->P()
	       <<"\n  after  "<<p_ijt->P()+p_kt->P()
	       <<"\n  p_"<<p_i->Id().front()<<" = "<<p_i->P()
	       <<"\n  p_"<<p_j->Id().front()<<" = "<<p_j->P()
	       <<"\n  p_"<<p_k->Id().front()<<" = "<<p_k->P()
	       <<"\n  p_{"<<p_ijt->Id()<<"} -> "<<p_ijt->P()
	       <<"\n  p_{"<<p_kt->Id()<<"} -> "<<p_kt->P()
	       <<"\n}"<<std::endl;
  }
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): m_type = "<<m_type
		 <<" {\n  m_z = "<<m_z<<", m_y = "<<m_y
		 <<"\n  p_"<<p_i->Id().front()<<" = "<<p_i->P()
		 <<"\n  p_"<<p_j->Id().front()<<" = "<<p_j->P()
		 <<"\n  p_"<<p_k->Id().front()<<" = "<<p_k->P()
		 <<"\n  p_{"<<p_ijt->Id()<<"} -> "<<p_ijt->P()
		 <<"\n  p_{"<<p_kt->Id()<<"} -> "<<p_kt->P()
		 <<"\n} -> "<<(m_type==2?1.0-m_y:m_y)
		 <<" vs. "<<p_info->AMin()<<" => stat = "
		 <<((m_type==2?1.0-m_y:m_y)>=p_info->AMin())<<std::endl;
#endif
}

void Dipole_Kinematics::EvaluateCSSKinematics()
{
  if (m_type==0) {
    double lrat=Lam(m_Q2,m_mij2,m_mk2)/Lam(m_Q2,(m_pi+m_pj).Abs2(),m_mk2);
    Vec4D pkt(sqrt(lrat)*(m_pk-(m_Q*m_pk/m_Q2)*m_Q)+
               (m_Q2+m_mk2-m_mij2)/(2.*m_Q2)*m_Q);
    p_ijt->SetP(m_Q-pkt);
    p_kt->SetP(pkt);
    double pijpk((m_pi+m_pj)*m_pk);
    m_z=(m_pi*m_pk)/pijpk;
    m_y=1.0/(1.0+pijpk/(m_pi*m_pj));
    if (Massive()) {
      double eps(m_Q2-m_mi2-m_mj2-m_mk2);
      m_ym=2.0*sqrt(m_mi2*m_mj2)/eps;
      m_yp=1.0-2.0*sqrt(m_mk2)*(sqrt(m_Q2)-sqrt(m_mk2))/eps;
      if (m_y<m_ym) p_info->SetStat(0);
    }
    PHASIC::Kin_Args ff(m_pi,m_pj,m_pk,m_y,m_z);
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,&ff,NULL);
    else m_kt2=(m_Q2-m_mi2-m_mj2-m_mk2)*m_y*m_z*(1.0-m_z)
              -sqr(1.0-m_z)*m_mi2-sqr(m_z)*m_mj2;
    if (p_info->Stat() && (m_pi[0]>1.0e-3 && m_pj[0]>1.0e-3) &&
        (pkt[0]<0.0 || m_Q[0]<pkt[0])) {
      p_info->SetStat(0);
      msg_Error()<<METHOD<<"(): Negative energy in FF {\n  p_i = "
                  <<m_pi<<"\n  p_j = "<<m_pj<<"\n  p_k = "<<m_pk
                  <<"\n  p_ij -> "<<m_Q-pkt<<"\n  p_k  -> "
                  <<pkt<<"\n}"<<std::endl;
    }
    for (size_t i(0);i<m_cur.size();++i) m_p[i]=m_cur[i]->P();
  }
  else if (m_type==2) {
    double pijpa((m_pi+m_pj)*m_pk);
    m_z=(m_pi*m_pk)/pijpa;
    m_y=-(m_pi*m_pj-0.5*(m_mij2-m_mi2-m_mj2))/pijpa;
    Vec4D pkt((1.0-m_y)*m_pk);
    p_ijt->SetP(m_Q-pkt);
    p_kt->SetP(pkt);
    if (Massive()) m_yp=1.0+m_y*(m_mij2-sqr(sqrt(m_mi2)+sqrt(m_mj2)))/m_Q2;
    PHASIC::Kin_Args ff(m_pi,m_pj,m_pk,1.0-m_y,m_z);
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,&ff,NULL);
    else m_kt2=2.0*(m_pi*m_pj)*m_z*(1.0-m_z)
              -sqr(1.0-m_z)*m_mi2-sqr(m_z)*m_mj2;
    if (p_info->Stat() && m_Q[0]<pkt[0] && m_pi[0]>1.0e-3 && m_pj[0]>1.0e-3) {
      p_info->SetStat(0);
      msg_Error()<<METHOD<<"(): Negative energy in FI {\n  p_i = "
                  <<m_pi<<"\n  p_j = "<<m_pj<<"\n  p_a = "<<m_pk
                  <<"\n  p_ij -> "<<m_Q-pkt<<"\n  p_a  -> "
                  <<pkt<<"\n}"<<std::endl;
    }
    for (size_t i(0);i<m_cur.size();++i) m_p[i]=m_cur[i]->P();
  }
  else if (m_type==1) {
    double pjpa=m_pj*m_pi, pkpa=m_pk*m_pi, pjpk=m_pj*m_pk;
    m_z=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
    m_y=pjpa/(pjpa+pkpa);
    p_ijt->SetP(m_z*m_pi);
    p_kt->SetP(m_Q-p_ijt->P());
    PHASIC::Kin_Args ff(m_pi,m_pj,m_pk,m_y,m_z);
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,&ff,NULL);
    else m_kt2=(-m_Q2+m_mk2)*m_y/m_z*(1.0-m_y)*(1.0-m_z);
    for (size_t i(0);i<m_cur.size();++i) m_p[i]=m_cur[i]->P();
  }
  else if (m_type==3) {
    double papb=m_pi*m_pk, pjpa=m_pj*m_pi, pjpb=m_pj*m_pk;
    m_z=(papb+pjpa+pjpb)/papb;
    Vec4D pajt(m_z*m_pi), K(-m_pi-m_pk-m_pj), Kt(-pajt-m_pk), KpKt(K+Kt);
    pajt=pajt-2.0*pajt*KpKt/(KpKt*KpKt)*KpKt+2.0*pajt*Kt/(K*K)*K;
    p_ijt->SetP(pajt);
    p_kt->SetP(m_pi+m_pj+m_pk-pajt);
    m_y=-pjpa/papb;
    PHASIC::Kin_Args ff(m_pi,m_pj,m_pk,m_y,m_z);
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,&ff,NULL);
    else m_kt2=m_Q2*m_y/m_z*(1.0-m_z-m_y);
    for (size_t i(0);i<m_cur.size();++i) {
      const Vec4D &p(m_cur[i]->P());
      m_p[i]=p-2.0*p*KpKt/(KpKt*KpKt)*KpKt+2.0*p*K/(Kt*Kt)*Kt;
    }
    m_pi=m_pi-2.0*m_pi*KpKt/(KpKt*KpKt)*KpKt+2.0*m_pi*Kt/(K*K)*K;
    m_pj=m_pj-2.0*m_pj*KpKt/(KpKt*KpKt)*KpKt+2.0*m_pj*Kt/(K*K)*K;
    m_pk=m_pk-2.0*m_pk*KpKt/(KpKt*KpKt)*KpKt+2.0*m_pk*Kt/(K*K)*K;
  }
  else {
    THROW(fatal_error,"Invalid dipole type");
  }
}

void Dipole_Kinematics::EvaluateAlaricKinematics()
{
  int index_i(p_i->Id().front()), index_j(p_j->Id().front()),
          index_k(p_k->Id().front());
  Cluster_Amplitude* ampl = Cluster_Amplitude::New();
  ampl->SetNIn(2);
  for (int i = 0; i <= m_cur.size(); ++i) {
    ampl->CreateLeg((i<2?-1.:1.)*p_subevt->p_real->p_mom[i],i<2?p_subevt->p_real->p_fl[i].Bar():p_subevt->p_real->p_fl[i]);
  }
  if (m_cur.back()->SubType()>>2&1) { // soft
    if (p_j->Flav().IsFermion()) std::swap(index_i, index_j);
    PHASIC::Ant_Args ff;
    for (int i = 0; i <= m_cur.size(); ++i) {
      ff.m_p.emplace_back((i<2?-1.:1.)*p_subevt->p_real->p_mom[i]);
    }

    ff.m_b=p_softrecoil->RecoilTags(ampl,(1<<index_i)|(1<<index_j),1<<index_k);
    PHASIC::ClusterAntenna(ff, index_i, index_j, index_k, 0.);

    p_ijt->SetP(ff.m_pijt);
    p_kt->SetP(ff.m_p[index_k]);
    m_n = ff.m_n;
    m_z = ff.m_z;
    m_y = ff.m_y;
  } else if (m_cur.back()->SubType()>>3&1) { // collinear
    Vec4D K = p_collrecoil->Recoil(ampl,(1<<index_i)|(1<<index_j),1<<index_k);
    std::vector<int> tags = p_collrecoil->RecoilTags(ampl,(1<<index_i)|(1<<index_j),1<<index_k);
    int nk = std::count_if(tags.begin(),tags.end(),[](int t){return t&2;});

    double K2(K.Abs2());
    int mode = 0;
    PHASIC::Kin_Args ff=PHASIC::ClusterFFDipole(0,0,0,K2,ampl->Mom(index_i),ampl->Mom(index_j),K,mode);
    if (ff.m_stat<0) {
      msg_Error()<<METHOD<<": Clustering failed in subtraction.\n";
    }

    if(nk>1) {
      Poincare oldcms(K), newcms(ff.m_pk);
      newcms.Invert();
      for(size_t i(0);i<ampl->Legs().size();++i) {
        if(tags[i]&2) {
          ampl->Leg(i)->SetMom(newcms*(oldcms*ampl->Leg(i)->Mom()));
        }
      }
    }
    else {
      for(size_t i(0);i<ampl->Legs().size();++i) {
        if(tags[i]&2) {
          ampl->Leg(i)->SetMom(ff.m_pk);
        }
      }
    }

    p_ijt->SetP(ff.m_pi);
    p_kt->SetP(ampl->Leg(index_k)->Mom());
    m_n=ff.m_nb;
    m_z=ff.m_z;
    m_y=ff.m_y;
  }
  else
    msg_Error() << METHOD << ": Unknown subtraction type, neither soft nor collinear.\n";
  if (m_type==0) {
    m_z=m_pi*m_n/((m_pi+m_pj)*m_n);
    m_y=m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);
    if (Massive()) {
      // TODO later
    }
    PHASIC::Kin_Args ff(m_pi,m_pj,m_pk,m_y,m_z);
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,&ff,NULL);
    else m_kt2=0.;
    if (p_info->Stat() && (m_pi[0]>1.0e-3 && m_pj[0]>1.0e-3) &&
        (p_kt->P()[0]<0.0 || m_Q[0]<p_kt->P()[0])) {
      p_info->SetStat(0);
      msg_Error()<<METHOD<<"(): Negative energy in FF {\n  p_i = "
                  <<m_pi<<"\n  p_j = "<<m_pj<<"\n  p_k = "<<m_pk
                  <<"\n  p_ij -> "<<m_Q-p_kt->P()<<"\n  p_k  -> "
                  <<p_kt->P()<<"\n}"<<std::endl;
    }
  }
  int it(0);
  for (int i = 0; i < ampl->Momenta().size(); ++i) {
    if (i == index_j) continue;
    if (i == index_i) m_p[it]=p_ijt->P();
    else if (i == index_k) m_p[it]=p_k->P();
    else m_p[it]=(i<2?-1.:+1.)*ampl->Leg(i)->Mom();
    ++it;
  }
  ampl->Delete();
}

void Dipole_Kinematics::CheckKT2Min()
{
  static double kt2c[2];
  static bool didsetkt2c{ false };
  if (!didsetkt2c) {
    auto pss = Settings::GetMainSettings()["SHOWER"];
    kt2c[0] = pss["FS_PT2MIN"].Get<double>();
    kt2c[1] = pss["IS_PT2MIN"].Get<double>();
    didsetkt2c = true;
  }
  if (m_kt2<kt2c[m_type&1]) m_a=1.0;
}

std::ostream &METOOLS::operator<<
  (std::ostream &str,const Dipole_Kinematics &k)
{
  return str<<*k.JI()<<","<<*k.JJ()<<"<->"<<*k.JK()<<" "<<k.Type();
}
