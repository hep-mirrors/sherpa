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

std::set<ATOOLS::Particle*> SHERPA::Remnant_Base::s_last[2];

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
  p_partner(NULL),
  m_emin(0.0) {}

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
  s_last[0].clear();
  s_last[1].clear();
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
  ATOOLS::ms::error_code error=shift.Scale();
  if (error!=ATOOLS::ms::no_error) {
    ATOOLS::msg.Info()<<"Remnant_Base::AdjustKinematics(): "<<error<<".\n"
		      <<"   Retry using new remnant pair."<<std::endl;
    while (ChooseLast()) {
      if (AdjustKinematics()) return true;
    }
    return false;
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
  double E=parton->Momentum()[0];
  if (E<0.0 || E>m_ebeam) {
    ATOOLS::msg.Error()<<"Remnant_Base::Extract("<<parton<<"): "
		       <<"Constituent energy out of range E = "
		       <<E<<"."<<std::endl;
    return false;
  }
  m_erem-=parton->Momentum()[0]+MinimalEnergy(parton->Flav());
  if (m_erem<=0.0) {
    msg_Tracking()<<"Remnant_Base::Extract(..): No remaining energy for "
		  <<parton->Flav()<<", p = "<<parton->Momentum()<<" -> E_min = "
		  <<(parton->Momentum()[0]+MinimalEnergy(parton->Flav()))<<std::endl;
    return false;
  }
  if (parton->Momentum()[0]<=m_emin) {
    msg_Tracking()<<"Remnant_Base::Extract(..): Energy exceeds minimum for "
		  <<parton->Flav()<<", p = "<<parton->Momentum()<<" <- E_min = "
		  <<m_emin<<std::endl;
    return false;
  }
  return true;
}

bool Remnant_Base::FindLast(const short unsigned int side)
{
  ATOOLS::Blob *cur=p_beamblob;
  if (side==1) cur=p_partner->p_beamblob;
  for (int j=0;j<cur->NOutP();++j) {
    if (cur->OutParticle(j)->DecayBlob()==NULL &&
	s_last[side].find(cur->OutParticle(j))==s_last[side].end()) {
      s_last[side].insert(cur->OutParticle(j));
      p_last[side]=cur->OutParticle(j);
      return true;
    }
  }
  return false;
}

bool Remnant_Base::ChooseLast()
{
  if (FindLast(1)) return true;
  if (FindLast(0)) {
    s_last[1].clear();
    return FindLast(1);
  }
  return false;
}
