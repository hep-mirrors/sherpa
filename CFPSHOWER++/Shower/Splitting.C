#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CFPSHOWER {
  size_t Splitting::s_cnt=0;
}

ostream & CFPSHOWER::operator<<(std::ostream &s,const splitting_mode::code & mode) {
  s<<"mode = ";
  switch (int(mode)) {
  case 1:
    s<<"coll"; break;
  case 2:
    s<<"soft"; break;
  case 0:
  default:
    s<<"diff"; break;
  }
  return s;
}

ostream & CFPSHOWER::operator<<(ostream &s,const Splitting & split) {
  s<<"Splitting(splitter = "<<split.GetSplitter()->Id()<<", "
   <<"spectator = "<<split.GetSpectator()->Id()<<", "
   <<"t = "<<split.t()<<", t2 = "<<split.t2()<<", "
   <<"z = "<<split.z()<<", z2 = "<<split.z2()<<", "
   <<split.Mode()<<", "
   <<(*split.GetKernel());
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & tcut) :
  p_splitter(splitter), p_spectator(spectator), p_kernel(NULL), m_mode(splitting_mode::diff),
  m_Q2(0.), m_Q2red(m_Q2), 
  m_tstart(t), m_tcut(tcut),
  m_t(t), m_t2(-1.), m_z(-1.), m_z2(-1.), m_phi(-1.), m_phi2(-1.), m_y(0.), 
  m_isend(false), 
  p_weight(NULL),
  m_kinscheme(-1), m_clustered(0)
{
  if (p_splitter!=NULL && p_spectator!=NULL) {
    m_Q2  = ATOOLS::dabs(((p_splitter->Beam()>0?1.:-1.)*p_splitter->Mom() +
			  (p_spectator->Beam()>0?1.:-1.)*p_spectator->Mom()).Abs2());
    m_eta = splitter->XB();
  }
  s_cnt++;
  //msg_Out()<<"   --- "<<METHOD<<" for |"<<this<<"| ---------------------------------------\n";
}

Splitting::~Splitting() {
  //msg_Out()<<"   --- "<<METHOD<<" for |"<<this<<"| ---------------------------------------\n";
  if (p_weight) { delete p_weight; p_weight = NULL; }
  s_cnt--;
}

void Splitting::InitSplitting(const ATOOLS::Mass_Selector * ms) {
  m_eta       = p_splitter->XB(); 
  m_msplit2   = sqr(ms->Mass(p_splitter->Flav()));
  m_mspect2   = sqr(ms->Mass(p_spectator->Flav()));
  m_Q2red     = m_Q2 - m_mspect2;
  m_ismassive = (m_mspect2>0.);
  m_outs.clear(); m_m2s.clear();
  if (p_kernel==0 || p_kernel->GetSF()==0) {
    msg_Error()<<METHOD<<" throws error: "<<p_kernel<<(p_kernel==0?0:p_kernel->GetSF())<<"\n";
    return;
  }
  m_nflavs = p_kernel->GetSF()->GetFlavs().size();
  m_outs.resize(m_nflavs,NULL); m_m2s.resize(m_nflavs, 0.);
  for (size_t i=0;i<m_nflavs;i++) {
    m_m2s[i]   = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[i]));
    m_Q2red   -= m_m2s[i];
    if (m_m2s[i]>0.) m_ismassive = true;
  }
}

void Splitting::UpdateSpectatorMomentum() {
  p_spectator->SetMom(p_kernel->GetSF()->SpecMom());
  // update Bjorken-x variables of particles, where appropriate
  if (p_kernel->GetType()==kernel_type::IF ||
      p_kernel->GetType()==kernel_type::II) p_splitter->SetXB();
  if (p_kernel->GetType()==kernel_type::FI ||
      p_kernel->GetType()==kernel_type::II) p_spectator->SetXB();
}


