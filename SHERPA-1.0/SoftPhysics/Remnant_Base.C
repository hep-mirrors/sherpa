#include "Remnant_Base.H"

#include "MathTools.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace SHERPA;

Remnant_Base::Remnant_Base(const TypeID _m_type,const unsigned int _m_beam):
  m_type(_m_type),
  m_beam(_m_beam),
  p_partner(NULL) {}

Remnant_Base::~Remnant_Base() {}

void Remnant_Base::Clear()
{
  for (size_t i=0;i<3;++i) m_parton[i].clear();
  m_active=true;
  p_last[1]=p_last[0]=NULL;
  p_beamblob=NULL;
}

bool Remnant_Base::AdjustKinematics()
{
  if (!m_active) return true;
  if (p_partner==NULL) {
    ATOOLS::msg.Error()<<"Remnant_Base::AdjustKinematics(): "
		       <<"No partner found. Abort."<<std::endl;
    exit(129);
  }
  p_last[1]=p_partner->Last();
  if ((p_last[0]==NULL)||(p_last[1]==NULL)) {
    ATOOLS::msg.Error()<<"Remnant_Base::AdjustKinematics(): "
		       <<"Not enough remnants to ensure four momentum conservation. Abort."<<std::endl;
    exit(129);
  }
  double Erem=ATOOLS::rpa.gen.Ecms(), pzrem=0.0;
  for (int i=0;i<p_beamblob->NOutP();++i) {
    if (p_beamblob->OutParticle(i)!=p_last[0]) {
      Erem-=p_beamblob->OutParticle(i)->Momentum()[0];
      pzrem-=p_beamblob->OutParticle(i)->Momentum()[3];
    }
  }
  for (int i=0;i<p_partner->BeamBlob()->NOutP();++i) {
    if (p_partner->BeamBlob()->OutParticle(i)!=p_last[1]) {
      Erem-=p_partner->BeamBlob()->OutParticle(i)->Momentum()[0];
      pzrem-=p_partner->BeamBlob()->OutParticle(i)->Momentum()[3];
    }
  }
  ATOOLS::Vec4D pr1=p_last[0]->Momentum(), pr2=p_last[1]->Momentum();
  double sp, sp1, sp2, lam2, c1, c2, c3, E1, E2, pz1, pz2;
  sp1=ATOOLS::sqr(p_last[0]->Flav().PSMass())+ATOOLS::sqr(pr1[1])+ATOOLS::sqr(pr1[2]);
  sp2=ATOOLS::sqr(p_last[1]->Flav().PSMass())+ATOOLS::sqr(pr2[1])+ATOOLS::sqr(pr2[2]);
  sp=Erem*Erem-pzrem*pzrem;
  lam2=(sp-sp1-sp2)*(sp-sp1-sp2)/4.0-sp1*sp2;
  c1=0.5*(sp-sp1+sp2); c2=0.5*(sp+sp1-sp2); 
  c3=0.5*(sp-sp1-sp2)*Erem*Erem-lam2;
  double spn, ytn, yto=(Erem+pzrem)/(Erem-pzrem);
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=Erem*c2/(c1+c2)*(1.0+sign*sqrt(1.0+(c3/(Erem*Erem)-c2)/(c2*c2)*(c1+c2)));
    E2=Erem-E1;
    pz1=ATOOLS::Sign(pr1[3])*sqrt(E1*E1-sp1); pz2=ATOOLS::Sign(pr2[3])*sqrt(E2*E2-sp2);
    spn=ATOOLS::sqr(Erem)-ATOOLS::sqr(pz1+pz2);
    if (ATOOLS::dabs((spn-sp)/(spn+sp))>ATOOLS::rpa.gen.Accu()) 
      { pz1*=-1.0; spn=ATOOLS::sqr(E1+E2)-ATOOLS::sqr(pz1+pz2); }
    if (ATOOLS::dabs((spn-sp)/(spn+sp))>ATOOLS::rpa.gen.Accu()) continue;
    ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2);
    if (ATOOLS::dabs((ytn-yto)/(ytn+yto))>ATOOLS::rpa.gen.Accu()) 
      { pz1*=-1.0; pz2*=-1.0; ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2); }
    if (ATOOLS::dabs((ytn-yto)/(ytn+yto))>ATOOLS::rpa.gen.Accu()) continue;
  }
  pr1=ATOOLS::Vec4D(E1,pr1[1],pr1[2],pz1); p_last[0]->SetMomentum(pr1);
  pr2=ATOOLS::Vec4D(E2,pr2[1],pr2[2],pz2); p_last[1]->SetMomentum(pr2);
  return true;
}

const Remnant_Base::TypeID Remnant_Base::Type() const
{ return m_type; }

const unsigned int Remnant_Base::Beam() const
{ return m_beam; }

ATOOLS::Particle *Remnant_Base::Last()
{ m_active=false; return p_last[0]; }

ATOOLS::Blob *Remnant_Base::BeamBlob()
{ return p_beamblob; }

void Remnant_Base::ExtractParton(ATOOLS::Particle *_m_parton)
{ if (_m_parton!=NULL) m_parton[1].push_back(_m_parton); }

void Remnant_Base::SetPartner(SHERPA::Remnant_Base *_p_partner)
{ p_partner=_p_partner; }

