#include "METOOLS/Explicit/Dipole_Kinematics.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;

Dipole_Kinematics::Dipole_Kinematics
(Dipole_Info *const info,Current *const i,Current *const j,
 Current *const k,Current *const ijt,Current *const kt):
  p_i(i), p_j(j), p_k(k), p_ijt(ijt), p_kt(kt),
  m_type(0), m_swap(0), m_trig(1), p_info(info), m_mi2(0.0), m_mj2(0.0),
  m_mij2(sqr(p_ijt->Flav().Mass())), m_mk2(sqr(p_k->Flav().Mass())),
  m_ym(0.0), m_yp(1.0), m_f(0.0), m_a(0.0)
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
  m_evol=ToType<int>(rpa->gen.Variable("CSS_EVOLUTION_SCHEME"));
  m_kt2c[0]=ToType<double>(rpa->gen.Variable("CSS_FS_PT2MIN"));
  m_kt2c[1]=ToType<double>(rpa->gen.Variable("CSS_IS_PT2MIN"));
}

Dipole_Kinematics::~Dipole_Kinematics()
{
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
  if (m_type==0) {
    double lrat=Lam(m_Q2,m_mij2,m_mk2)/Lam(m_Q2,(m_pi+m_pj).Abs2(),m_mk2);
    Vec4D pkt(sqrt(lrat)*(m_pk-(m_Q*m_pk/m_Q2)*m_Q)+
	      (m_Q2+m_mk2-m_mij2)/(2.*m_Q2)*m_Q);
    if (Massive()) {
      double eps(m_Q2-m_mi2-m_mj2-m_mk2);
      m_ym=2.0*sqrt(m_mi2*m_mj2)/eps;
      m_yp=1.0-2.0*sqrt(m_mk2)*(sqrt(m_Q2)-sqrt(m_mk2))/eps;
    }
    p_ijt->SetP(m_Q-pkt);
    p_kt->SetP(pkt);
    double pijpk((m_pi+m_pj)*m_pk);
    m_z=(m_pi*m_pk)/pijpk;
    m_y=1.0/(1.0+pijpk/(m_pi*m_pj));
    if (m_evol==0) {
      m_kt2=(m_Q2-m_mi2-m_mj2-m_mk2)*m_y*m_z*(1.0-m_z)-
	sqr(1.0-m_z)*m_mi2-sqr(m_z)*m_mj2;
    }
    else {
      m_kt2=2.0*(m_pi*m_pj)*m_z*(1.0-m_z);
      if (p_ijt->Flav().IsFermion())
	m_kt2=2.0*(m_pi*m_pj)*(p_i->Flav().IsFermion()?(1.0-m_z):m_z);
      else if (p_i->Flav().IsFermion()) m_kt2=2.0*(m_pi*m_pj);
    }
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
    double lrat=Lam(m_Q2,m_mij2,m_mk2)/Lam(m_Q2,(m_pi+m_pj).Abs2(),m_mk2);
    Vec4D pkt(sqrt(lrat)*(m_pk-(m_Q*m_pk/m_Q2)*m_Q)+
	      (m_Q2+m_mk2-m_mij2)/(2.*m_Q2)*m_Q);
    p_ijt->SetP(m_Q-pkt);
    p_kt->SetP(pkt);
    double pijpa((m_pi+m_pj)*m_pk);
    m_z=(m_pi*m_pk)/pijpa;
    m_y=-(m_pi*m_pj)/pijpa;
    if (Massive()) m_yp=1.0+m_y*(m_mij2-sqr(sqrt(m_mi2)+sqrt(m_mj2)))/m_Q2;
    if (m_evol==0) {
      m_kt2=2.0*(m_pi*m_pj)*m_z*(1.0-m_z)-sqr(1.0-m_z)*m_mi2-sqr(m_z)*m_mj2;
    }
    else {
      m_kt2=2.0*(m_pi*m_pj)*m_z*(1.0-m_z);
      if (p_ijt->Flav().IsFermion())
	m_kt2=2.0*(m_pi*m_pj)*(p_i->Flav().IsFermion()?(1.0-m_z):m_z);
      else if (p_i->Flav().IsFermion()) m_kt2=2.0*(m_pi*m_pj);
    }
    if (p_info->Stat() && (m_pi[0]>1.0e-3 && m_pj[0]>1.0e-3) &&
	(m_Q[0]<pkt[0] || lrat<0.0)) {
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
    double sjk=(m_pj+m_pk).Abs2();
    double lrat=Lam(m_Q2,m_mk2,0.0)/Lam(m_Q2,sjk,0.0);
    p_ijt->SetP(sqrt(lrat)*(m_pi-(m_Q2-sjk)/(2.0*m_Q2)*m_Q)
		+(m_Q2-m_mk2)/(2.0*m_Q2)*m_Q);
    p_kt->SetP(m_Q-p_ijt->P());
    m_z=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
    m_y=pjpa/(pjpa+pkpa);
    if (m_evol==0) {
      m_kt2=(-m_Q2+m_mk2)*m_y/m_z*(1.0-m_z);
    }
    else {
      m_kt2=(-m_Q2+m_mk2)*m_y/m_z*(1.0-m_z);
      if (p_j->Flav().IsFermion()) m_kt2=(-m_Q2+m_mk2)*m_y/m_z;
    }
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
    if (m_evol==0) {
      m_kt2=m_Q2*m_y/m_z*(1.0-m_z-m_y);
    }
    else {
      m_kt2=m_Q2*m_y/m_z*(1.0-m_z-m_y);
      if (p_j->Flav().IsFermion()) m_kt2=m_Q2*m_y/m_z;
    }
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
  m_trig=m_y/m_yp<p_info->AMax(m_type);
  if (m_trig) m_trig=m_kt2<p_info->KT2Max();
  if (m_y<Max(m_ym,p_info->AMin())) p_info->SetStat(0);
  if (p_info->Stat() &&
      p_i->P()+p_j->P()+p_k->P()!=p_ijt->P()+p_kt->P()) {
    msg_Error()<<METHOD<<"(): Momentum not conserved in type "<<m_type
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

void Dipole_Kinematics::CheckKT2Min()
{
  if (m_kt2<m_kt2c[m_type&1]) m_a=1.0;
}

std::ostream &METOOLS::operator<<
  (std::ostream &str,const Dipole_Kinematics &k)
{
  return str<<*k.JI()<<","<<*k.JJ()<<"<->"<<*k.JK()<<" "<<k.Type();
}
