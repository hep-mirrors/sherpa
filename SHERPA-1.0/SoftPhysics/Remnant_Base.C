#include "Remnant_Base.H"

#include "Run_Parameter.H"
#include "Exception.H"
#include "Momentum_Shifter.H"
#include "MI_Handler.H"
#include "Scaling.H"

#ifdef PROFILE__all
#define PROFILE__Remnant_Base
#endif
#ifdef PROFILE__Remnant_Base
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

std::ostream &SHERPA::operator<<(std::ostream &ostr,const rtp::code code)
{
  switch (code) {
  case rtp::intact:      return ostr<<"Intact";
  case rtp::qcd_remnant: return ostr<<"QCD Remnant";
  case rtp::hadron:      return ostr<<"Hadron";
  case rtp::photon:      return ostr<<"Photon";
  case rtp::electron:    return ostr<<"Electron";
  }
  return ostr;
}

Remnant_Base::Remnant_Base(const rtp::code type,const unsigned int beam):
  Object("Remnant_Base_"+ATOOLS::ToString(beam)),
  m_type(type),
  m_beam(beam),
  p_partner(NULL) {}

Remnant_Base::~Remnant_Base() {}

void Remnant_Base::Clear()
{
  m_extracted.clear();
  m_companions.clear();
  m_active=true;
  p_last[1]=p_last[0]=NULL;
  p_beamblob=NULL;
  m_erem=m_ebeam;
  m_initialized=false;
}

void Remnant_Base::QuickClear()
{
  m_extracted.clear();
  m_erem=m_ebeam;
  m_initialized=false;
}

bool Remnant_Base::AdjustKinematics()
{
  PROFILE_HERE;
  if (!m_active) return true;
  if (p_partner==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,
			    "No partner remnant found.",
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
    shift.SetShift(ATOOLS::Vec4D(m_erem-(pr1+pr2)[0],0.,0.,
				 m_pzrem-(pr1+pr2)[3]));
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
			 <<"Parton ("<<p_last[i]<<") has non-positive energy "
			 <<p_last[i]->Momentum()<<std::endl;
      return false;
    }
  }
  return true;
}

void Remnant_Base::UnDo() 
{
  msg_Tracking()<<"Remnant_Base::UnDo(): "
		<<"Undoing changes on blob list."<<std::endl;
  while (p_beamblob->NOutP()>0) {
    p_beamblob->RemoveOutParticle(p_beamblob->OutParticle(0));
  }
  while (m_companions.size()>0) {
    delete *m_companions.begin();
    m_companions.erase(m_companions.begin());
  }
  ++m_errors;
}

double Remnant_Base::MinimalEnergy(const ATOOLS::Flavour &flavour) 
{
  return 0.;
}

ATOOLS::Flavour Remnant_Base::ConstituentType(const ATOOLS::Flavour &flavour) 
{
  return ATOOLS::kf::none;
}

bool Remnant_Base::Extract(ATOOLS::Particle *parton) 
{ 
  m_extracted.push_back(parton); 
  m_erem-=parton->Momentum()[0]+MinimalEnergy(parton->Flav());
  if (m_erem<=0.0) {
    msg_Tracking()<<"Remnant_Base::Extract(..): No remaining energy for "<<parton->Flav()
		  <<", p = "<<parton->Momentum()<<" -> E_min = "
		  <<(parton->Momentum()[0]+MinimalEnergy(parton->Flav()))<<std::endl;
    return false;
  }
  return true;
}

bool Remnant_Base::FindHardProcess(ATOOLS::Particle *const initiator,
				   ATOOLS::Blob *&process)
{
  ATOOLS::Blob *cur=initiator->DecayBlob();
  if (cur!=NULL) {
    if (cur->Type()==ATOOLS::btp::Hard_Collision) {
      process=cur;
      return true;
    }
    else {
      for (size_t i=0;i<(size_t)cur->NOutP();++i) {
	if (FindHardProcess(cur->OutParticle(i),process)) return true;
      }
    }
  }
  return false;
}

void Remnant_Base::FindIncoming(ATOOLS::Blob *const blob,const size_t beam,
				ATOOLS::Particle *&incoming)
{
  PROFILE_HERE;
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

bool Remnant_Base::AcquireMass(const ATOOLS::Particle *left,
			       const ATOOLS::Particle *right,
			       const double newsp)
{ 
  PROFILE_HERE;
  ATOOLS::Vec4D P=left->Momentum()+right->Momentum();
  ATOOLS::Vec4D p=p_last[0]->Momentum()+p_last[1]->Momentum();
  double C=(P[0]+p[0])/(P[3]+p[3]);
  double D=0.5*(p.MPerp2()-newsp)/(P[3]+p[3]);
  double E=P[0]-C*(P[3]-D);
  m_deltae=(E-sqrt(E*E+(1.0-C*C)*D*(D-2.0*P[3])))/(1.0-C*C);
  if (ATOOLS::Sign(P[3]-m_deltap)!=ATOOLS::Sign(P[3])) m_deltap*=-1.0;
  if (!(m_deltae>0.0) || P[0]<m_deltae) return false;
  // if (p_mihandler!=NULL) if (m_deltae>p_mihandler->ScaleMin(0)) return false;
  m_deltap=m_deltae*C+D;
  return true;
}

bool Remnant_Base::AdjustEnergy()
{
  PROFILE_HERE;
  ATOOLS::Blob *process=NULL;
  double sp1=ATOOLS::sqr(p_last[0]->Flav().PSMass())+
    p_last[0]->Momentum().PPerp2();
  double sp2=ATOOLS::sqr(p_last[1]->Flav().PSMass())+
    p_last[1]->Momentum().PPerp2();
  double spmin=2.*sqrt(sp1*sp2)+sp1+sp2;
  for (size_t i=0;i<m_extracted.size();++i) {
    if (!FindHardProcess(m_extracted[i],process)) return false;
    ATOOLS::Particle *in[2];
    for (short unsigned int i=0;i<2;++i) FindIncoming(process,i,in[i]);
    if (AcquireMass(in[0],in[1],spmin)) {
      ATOOLS::Momentum_Shifter boost(in[0],in[1]);
      boost.SetShift(ATOOLS::Vec4D(-m_deltae,0.,0.,-m_deltap));
      if (boost.Boost()!=ATOOLS::ms::no_error) continue;
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

