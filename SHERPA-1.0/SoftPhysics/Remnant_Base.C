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
  shift.SetShift(ATOOLS::Vec4D(m_erem-(pr1+pr2)[0],0.,0.,m_pzrem-(pr1+pr2)[3]));
  if (shift.Scale()!=ATOOLS::ms::no_error) AdjustEnergy();
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

void Remnant_Base::AdjustEnergy()
{
}
