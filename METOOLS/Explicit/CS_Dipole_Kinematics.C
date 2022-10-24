#include "METOOLS/Explicit/CS_Dipole_Kinematics.H"

#include "METOOLS/Explicit/Vertex.H"
#include "PDF/Main/NLOMC_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace METOOLS;
using namespace ATOOLS;

void CS_Dipole_Kinematics::Evaluate()
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
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,m_z,m_y,m_Q2);
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
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,m_z,1.0-m_y,m_Q2);
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
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,m_z,m_y,m_Q2);
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
    if (p_nlomc) m_kt2=p_nlomc->KT2(*p_subevt,m_z,m_y,m_Q2);
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
  m_trig = m_y/m_yp < ((p_nlomc) ? 1.0 : p_info->AMax(m_type));
  if (m_trig) m_trig=m_kt2<p_info->KT2Max();
  if (p_info->AMin()>0.0) {
    if (m_y<p_info->AMin()) p_info->SetStat(0);
  }
  else {
    if (m_kt2<-p_info->AMin()) p_info->SetStat(0);
  }
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


void CS_Dipole_Kinematics::EvaluateFFS(const Vec4D &p,Vec4D &q,
                                    double &Ai,double &B,double &t) {
  Ai = 0;
  B = 0;
  t = 0;
  q = Vec4D();
  if (Type()==0) {
    double zi(Z()), zj(1.0-zi), rv(1.0);
    Vec4D pi(PI()), pj(PJ());
    double pij2((pi+pj).Abs2()), mt2(0.0);
    double zim(zi), zjm(zj), zti(zi), ztj(zj);
    if (p_info->SubType()==1) {
      zti=1.0-(1.0-zti)*(1.0-Y());
      ztj=1.0-(1.0-ztj)*(1.0-Y());
    }
    SetA(0.5*(zti*(1.0-zti)+ztj*(1.0-ztj)));
    if (Massive()) {
      double y(Y()), q2(Q2());
      double s(q2-m_mi2-m_mj2-m_mk2);
      double viji(sqrt(sqr(s*y)-sqr(2.0*m_mi*m_mj))/(s*y+2.0*m_mi2));
      rv=sqrt(sqr(2.0*m_mk2+s*(1.0-y))-4.0*m_mk2*q2)/(s*(1.0-y));
      double zc(0.5*(2.0*m_mi2+s*y)/(m_mi2+m_mj2+s*y));
      double zm(zc*(1.0-viji*rv)), zp(zc*(1.0+viji*rv));
      mt2=2.0*p_info->Kappa()*(zp*zm-m_mi2/pij2);
      zim-=0.5*(1.0-rv);
      zjm-=0.5*(1.0-rv);
      SetA(A()-zp*zm);
    }
    q=zim*pi-zjm*pj;
    Ai=1.0-mt2;
    B=4.0*A();
    t=rv*pij2;
    SetA(Ai-2.*A());
  }
  else if (Type()==2) {
    double zi(Z()), zj(1.0-zi);
    Vec4D pi(PI()), pj(PJ());
    double pij2((pi+pj).Abs2());
    q=zi*pi-zj*pj;
    Ai=1.0;
    B=-4.0*q.Abs2()/pij2;
    t=pij2*(1.0-Y());
    SetA((1.0-zi)*zi);
    if (Massive()) {
      double Q2((pi+pj+PK()).Abs2());
      double mui2(m_mi2/Q2), y(Y());
      double eps(sqrt(sqr(y-2.0*mui2)-4.0*mui2*mui2)/y);
      SetA((0.5*(1.0+eps)-zi)*(zi-0.5*(1.0-eps)));
    }
    SetA(Ai-2.0*A());
  }
  else if (Type()==1) {
    double x(Z()), ui(Y());
    Vec4D pi(PJ()), pk(PK());
    q=pi/ui-pk/(1.0-ui);
    Ai=x;
    double tc((1.0-x)/x);
    B=2.0*tc*ui*(1.0-ui)*q.Abs2()/(pi*pk);
    t=-2.0*(pi*PI())*x;
    SetA(tc);
    if (Massive()) {
      double Q2(2.0*(JKT()->P()*p));
      SetA(tc-pk.Abs2()/Q2*ui/(1.0-ui));
    }
    SetA(Ai+2.0*A());
  }
  else {
    double x(Z()), vi(Y());
    Vec4D pi(PJ()), pk(-PK());
    double z(x), tc((1.0-x)/x);
    if (p_info->SubType()&3) z=x+vi;
    if (p_info->SubType()&3) tc+=1.0/(x+vi)-1.0/x;
    Ai=z;
    B=-4.0*tc;
    q=pi-vi*pk;
    t=-2.0*(pi*PI())*x;
    SetA(Ai+2.0*tc);
  }
}

void CS_Dipole_Kinematics::EvaluateFVS(double &Ai,double &t) {
  Ai = 0;
  t = 0;
  bool iisf(JI()->Flav().IsFermion());
  if (Type()==0) {
    double zi(Z()), y(Y());
    double zti(zi), ztj(1.0-zi);
    if (p_info->SubType()==1) {
      zti=1.0-(1.0-zti)*(1.0-y);
      ztj=1.0-(1.0-ztj)*(1.0-y);
    }
    double pipj(PI()*PJ());
    double rv(1.0), mt2(0.0);
    if (Massive()) {
      double pij2(2.0*pipj+m_mij2), q2(Q2());
      rv=(q2-pij2-m_mk2)/(q2-m_mij2-m_mk2)*
        sqrt((sqr(q2-m_mij2-m_mk2)-4.0*m_mij2*m_mk2)/
             (sqr(q2-pij2-m_mk2)-4.0*pij2*m_mk2));
      mt2=m_mij2/pipj;
    }
    if (iisf) Ai=2.0/(1.0-zi*(1.0-y))-rv*(1.0+zti+mt2);
    else Ai=2.0/(1.0-(1.0-zi)*(1.0-y))-rv*(1.0+ztj+mt2);
    t=2.0*pipj;
    SetA(Ai);
  }
  else if (Type()==2) {
    double zi(Z()), y(Y());
    double pipj(PI()*PJ());
    double mt2(m_mij?m_mij2/pipj:0.0);
    if (iisf) Ai=2.0/(1.0-zi+y)-(1.0+zi+mt2);
    else Ai=2.0/(1.0-(1.0-zi)+y)-(2.0-zi+mt2);
    if (p_info->SubType()==2 &&
        !Massive()) {
      if (iisf) Ai=2.0*zi/(1.0-zi+y)+(1.0-zi);
      else Ai=2.0*(1.0-zi)/(1.0-(1.0-zi)+y)+zi;
    }
    t=2.0*pipj*(1.0-y);
    SetA(Ai);
  }
  else if (Type()==1) {
    double x(Z()), ui(Y());
    if (iisf) Ai=2.0/(1.0-x+ui)-(1.0+x);
    else Ai=1.0-2.0*x*(1.0-x);
    if (p_info->SubType()==2 &&
        !Massive()) {
      if (iisf) Ai=2.0*x/(1.0-x+ui)+(1.0-x);
    }
    t=-2.0*(PI()*PJ())*x;
    SetA(Ai);
  }
  else {
    double x(Z()), z(x);
    if (p_info->SubType()&3) z=x+Y();
    if (iisf) Ai=2.0/(1.0-x)-(1.0+z);
    else Ai=1.0-2.0*z*(1.0-z);
    if (p_info->SubType()==2) {
      if (iisf) Ai=2.0*z/(1.0-x)+(1.0-z);
      else Ai=1.0-2.0*z*(1.0-z);
    }
    t=-2.0*(PI()*PJ())*x;
    SetA(Ai);
  }
}


void CS_Dipole_Kinematics::EvaluateVVS(const ATOOLS::Vec4D &p, ATOOLS::Vec4D &q,
                                    double &Ai,double &Aj,double &B,double &t) {
  Ai = 0;
  Aj = 0;
  B = 0;
  t = 0;
  q = Vec4D();
  if (Type()==0) {
    double zi(Z()), zj(1.0-zi), y(Y());
    Vec4D pi(PI()), pj(PJ());
    double pipj(pi*pj), sl(1.0);
    double zim(zi), zjm(zj), zti(zi), ztj(zj);
    if (p_info->SubType()==1) {
      zti=1.0-(1.0-zti)*(1.0-Y());
      ztj=1.0-(1.0-ztj)*(1.0-Y());
    }
    SetA(0.5*(zti*(1.0-zti)+ztj*(1.0-ztj)));
    if (Massive()) {
      double y(Y()), q2(Q2()), s(q2-m_mk2);
      double rv(sqrt(sqr(2.0*m_mk2+s*(1.0-y))-4.0*m_mk2*q2)/(s*(1.0-y)));
      double zm(0.5*(1.0-rv)), zp(0.5*(1.0+rv));
      sl=(1.0-0.5*p_info->Kappa()*zp*zm)/rv;
      zim-=0.5*(1.0-rv);
      zjm-=0.5*(1.0-rv);
      SetA((A()-zp*zm)/rv);
    }
    Ai=2.0*(1.0/(1.0-zi*(1.0-y))-sl);
    Aj=2.0*(1.0/(1.0-zj*(1.0-y))-sl);
    if (Swap()) std::swap<double>(Ai,Aj);
    q=zim*pi-zjm*pj;
    B=-2.0*A();
    t=2.0*pipj;
    SetA(Ai+Aj+2.0*A());
  }
  else if (Type()==2) {
    double zi(Z()), zj(1.0-zi), y(Y());
    Vec4D pi(PI()), pj(PJ());
    Ai=2.0*(zi-y)/(1.0-zi+y);
    Aj=2.0*(zj-y)/(1.0-zj+y);
    if (p_info->SubType()==2) {
      Ai=2.0*zi/(1.0-zi+y);
      Aj=2.0*zj/(1.0-zj+y);
    }
    if (Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*zi*zj;
    q=zi*pi-zj*pj;
    t=2.0*(pi*pj)*(1.0-y);
    SetA(Ai+Aj+2.0*zi*zj);
  }
  else if (Type()==1) {
    double x(Z()), ui(Y());
    Vec4D pi(PJ()), pk(PK());
    Ai=2.0*(x-ui)/(1.0-x+ui);
    Aj=2.0*x*(1.0-x);
    if (p_info->SubType()==2 &&
	!Massive()) {
      Ai=2.0*x/(1.0-x+ui);
    }
    if (Swap()) std::swap<double>(Ai,Aj);
    q=pi/ui-pk/(1.0-ui);
    double tc((1.0-x)/x);
    B=tc*ui*(1.0-ui)*q.Abs2()/(pi*pk);
    t=-2.0*(pi*PI())*x;
    SetA(tc);
    if (Massive()) {
      double q2(2.0*(JKT()->P()*p));
      SetA(tc-pk.Abs2()/q2*ui/(1.0-ui));
    }
    SetA(Ai+Aj+2.0*A());
  }
  else {
    double x(Z()), vi(Y());
    Vec4D pi(PJ()), pk(-PK());
    double z(x), tc((1.0-x)/x);
    if (p_info->SubType()&3) z=x+vi;
    if (p_info->SubType()&3) tc+=1.0/(x+vi)-1.0/x;
    Ai=2.0*x/(1.0-x);
    Aj=2.0*z*(1.0-z);
    if (p_info->SubType()==2) Ai+=2.0*(z/(1.0-x)-x/(1.0-x));
    if (Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*tc;
    q=pi-vi*pk;
    t=-2.0*(pi*PI())*x;
    SetA(Ai+Aj+2.0*tc);
  }
}

void CS_Dipole_Kinematics::CheckKT2Min()
{
  static double kt2c[2];
  static bool didsetkt2c{ false };
  if (!didsetkt2c) {
    Settings& s = Settings::GetMainSettings();
    kt2c[0] = s["CSS_FS_PT2MIN"].Get<double>();
    kt2c[1] = s["CSS_IS_PT2MIN"].Get<double>();
    didsetkt2c = true;
  }
  if (m_kt2<kt2c[m_type&1]) m_a=1.0;
}


