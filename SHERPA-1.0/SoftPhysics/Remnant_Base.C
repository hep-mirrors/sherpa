#include "Remnant_Base.H"

#include "Run_Parameter.H"
#include "Exception.H"
#include "Momentum_Shifter.H"

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
  m_erem=ATOOLS::rpa.gen.Ecms();
  m_pzrem=0.0;
  for (size_t i=0;i<2;++i) {
    ATOOLS::Blob *cur=p_beamblob;
    if (i==1) cur=p_partner->p_beamblob;
    for (int j=0;j<cur->NOutP();++j) {
      if (cur->OutParticle(j)!=p_last[i]) {
	const ATOOLS::Vec4D &p=cur->OutParticle(j)->Momentum();
	m_erem-=p[0];
	m_pzrem-=p[3];
      }
    }
  }
  ATOOLS::Vec4D pr1=p_last[0]->Momentum(), pr2=p_last[1]->Momentum();
  ATOOLS::Momentum_Shifter shift(p_last[0],p_last[1]);
  shift.SetSPerp(ATOOLS::sqr(p_last[0]->Flav().PSMass())+pr1.PPerp2(),1);
  shift.SetSPerp(ATOOLS::sqr(p_last[1]->Flav().PSMass())+pr2.PPerp2(),2);
  if (!ATOOLS::IsZero(m_pzrem-(pr1+pr2)[3])) {
    shift.SetShift(ATOOLS::Vec4D(m_erem-(pr1+pr2)[0],0.,0.,m_pzrem-(pr1+pr2)[3]));
  }
  else {
    shift.SetDirection(ATOOLS::Vec4D(0.,0.,0.,1.));
  }
  if (shift.Scale()!=ATOOLS::ms::no_error) {
    if (!AdjustEnergy()) return false;
  }
  for (size_t i=0;i<2;++i) {
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

void Remnant_Base::FindHardProcess(ATOOLS::Particle *const initiator,ATOOLS::Blob *&process)
{
  ATOOLS::Blob *cur=initiator->DecayBlob();
  if (cur!=NULL) {
    if (cur->Type()==ATOOLS::btp::Signal_Process ||
	cur->Type()==ATOOLS::btp::Hard_Collision) {
      process=cur;
      return;
    }
    else {
      for (size_t i=0;i<(size_t)cur->NOutP();++i) FindHardProcess(cur->OutParticle(i),process);
    }
  }
}

void Remnant_Base::FindIncoming(ATOOLS::Blob *const blob,const size_t beam,
				ATOOLS::Particle *&incoming)
{
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    ATOOLS::Particle *cur=blob->InParticle(i);
    if (cur->ProductionBlob()->Type()==ATOOLS::btp::Beam) {
      if ((size_t)blob->Beam()==beam) incoming=cur;
    }
    else {
      FindIncoming(cur->ProductionBlob(),beam,incoming);
    }
  }
}

bool Remnant_Base::AcquireMass(ATOOLS::Blob *const process,const double newsp)
{
  ATOOLS::Vec4D P=process->InParticle(0)->Momentum()+process->InParticle(1)->Momentum();
  ATOOLS::Vec4D p=p_last[0]->Momentum()+p_last[1]->Momentum();
  double C=(P[0]+p[0])/(P[3]+p[3]), D=0.5*(p.MPerp2()-newsp)/(P[3]+p[3]), E=P[0]-C*(P[3]-D);
  m_deltae=(E-sqrt(E*E+(1.0-C*C)*D*(D-2.0*P[3])))/(1.0-C*C);
  if (!(m_deltae>0.0) || P[0]<m_deltae || 
      ATOOLS::Sign(P[3]-m_deltae)!=ATOOLS::Sign(P[3])) return false;
  m_deltap=m_deltae*C+D;
  return true;
}

bool Remnant_Base::AdjustEnergy()
{
  ATOOLS::Blob *process=NULL;
  double sp1=ATOOLS::sqr(p_last[0]->Flav().PSMass())+p_last[0]->Momentum().PPerp2();
  double sp2=ATOOLS::sqr(p_last[1]->Flav().PSMass())+p_last[1]->Momentum().PPerp2();
  double spmin=2.*sqrt(sp1*sp2)+sp1+sp2;
  for (size_t i=0;i<m_parton[1].size();++i) {
    FindHardProcess(m_parton[1][i],process);
    if (AcquireMass(process,spmin)) {
      ATOOLS::Particle *in[2];
      in[0]=process->InParticle(0);
      in[1]=process->InParticle(1);
      ATOOLS::Momentum_Shifter boost(in[0],in[1]);
      boost.SetShift(ATOOLS::Vec4D(-m_deltae,0.,0.,-m_deltap));
      if (boost.Boost()!=ATOOLS::ms::no_error) continue;
      FindIncoming(process,0,in[0]);
      FindIncoming(process,1,in[1]);
      boost.Boost(in[0]);
      boost.Boost(in[1]);
      ATOOLS::Vec4D pr1=p_last[0]->Momentum(), pr2=p_last[1]->Momentum();
      ATOOLS::Momentum_Shifter shift(p_last[0],p_last[1]);
      shift.SetSPerp(sp1,1);
      shift.SetSPerp(sp2,2);
      shift.SetShift(ATOOLS::Vec4D(m_erem+m_deltae-(pr1+pr2)[0],0.,0.,
				   m_pzrem+m_deltap-(pr1+pr2)[3]));
      if (shift.Scale()!=ATOOLS::ms::no_error) {
	boost.BoostBack();
	boost.BoostBack(in[0]);
	boost.BoostBack(in[1]);
	continue;
      }
      return true;
    }
  }
  return false;
}

