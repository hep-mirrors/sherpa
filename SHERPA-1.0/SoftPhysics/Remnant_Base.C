#include "Remnant_Base.H"

#include "MathTools.H"
#include "Run_Parameter.H"
#include "Exception.H"

using namespace SHERPA;

Remnant_Base::Remnant_Base(const TypeID type,const unsigned int beam):
  m_type(type),
  m_beam(beam),
  p_partner(NULL) {}

Remnant_Base::~Remnant_Base() {}

void Remnant_Base::Clear()
{
  for (size_t i=0;i<3;++i) m_parton[i].clear();
  m_active=true;
  p_last[1]=p_last[0]=NULL;
  p_beamblob=NULL;
  m_initialized=false;
}

double Remnant_Base::Lambda2(double sp,double sp1,double sp2) 
{ 
  return (sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
}

bool Remnant_Base::AdjustKinematics()
{
  if (!m_active) return true;
  if (p_partner==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"No partner remnant found.",
			    "Remnant_Base","AdjustKinematics"));
  }
  p_last[1]=p_partner->Last();
  if ((p_last[0]==NULL)||(p_last[1]==NULL)) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,
			    "Not enough remnants to ensure four momentum conservation.",
			    "Remnant_Base","AdjustKinematics"));
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
  double sp, sp1, sp2, E1, E2, pz1, pz2;
  sp1=ATOOLS::sqr(p_last[0]->Flav().PSMass())+ATOOLS::sqr(pr1[1])+ATOOLS::sqr(pr1[2]);
  sp2=ATOOLS::sqr(p_last[1]->Flav().PSMass())+ATOOLS::sqr(pr2[1])+ATOOLS::sqr(pr2[2]);
  sp=Erem*Erem-pzrem*pzrem;
  double spn, ytn, yto=(Erem+pzrem)/(Erem-pzrem);
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=0.5/sp*((sp+sp1-sp2)*Erem+sign*sqrt(Lambda2(sp,sp1,sp2))*pzrem);
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
    if (ATOOLS::dabs(pz1)>ATOOLS::dabs(pz2)) { if (ATOOLS::Sign(pz1)==ATOOLS::Sign(pr1[3])) break; }
    else { if (ATOOLS::Sign(pz2)==ATOOLS::Sign(pr2[3])) break; }
  }
  /*
  if (pz2>0.0&&pz1<0.0) {
    ATOOLS::msg.Error()<<"Remnant_Base::AdjustKinematics(): "
		       <<"Warning. Interchanged remnant momenta."<<std::endl
		       <<"   Former momenta are "<<pr1<<" "<<pr2<<std::endl
		       <<"   New momenta are    "<<ATOOLS::Vec4D(E1,pr1[1],pr1[2],pz1)
		       <<" "<<ATOOLS::Vec4D(E2,pr2[1],pr2[2],pz2)<<std::endl;
  }
  */
  pr1=ATOOLS::Vec4D(E1,pr1[1],pr1[2],pz1); p_last[0]->SetMomentum(pr1);
  pr2=ATOOLS::Vec4D(E2,pr2[1],pr2[2],pz2); p_last[1]->SetMomentum(pr2);
  for (size_t i=0;i<2;++i) {
    // the brackets are necessary for 'nan'-values
    if (!(p_last[i]->Momentum()[0]>0.)) {
      ATOOLS::msg.Error()<<"Remnant_Base::AdjustKinematics(): "
			 <<"Parton ("<<p_last[i]->Number()<<") has non-positive energy "
			 <<p_last[i]->Momentum()<<std::endl;
      UnDo();
      p_partner->UnDo();
      return false;
    }
  }
  return true;
}

void Remnant_Base::UnDo() 
{
  ATOOLS::msg.Tracking()<<"Remnant_Base::UnDo(): Undoing changes on blob list."<<std::endl;
  while (p_beamblob->NOutP()>0) {
    p_beamblob->RemoveOutParticle(p_beamblob->OutParticle(0));
  }
  for (ATOOLS::Particle_List::iterator pit=m_parton[0].begin();
       pit!=m_parton[0].end();++pit) {
    delete *pit;
  }
  m_parton[0].clear();
  m_parton[2].clear();
  ++m_errors;
}

double Remnant_Base::MinimalEnergy(const ATOOLS::Flavour &flavour) 
{
  return 0.;
}
