#include "CFPSHOWER++/Tools/Splitting.H"
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
   <<"t1 = "<<split.t(0)<<", t2 = "<<split.t(1)<<", "
   <<"z1 = "<<split.z(0)<<", z2 = "<<split.z(1)<<",\n"
   <<"    mode = "<<split.Mode()<<", "<<(*split.GetKernel());
  for (size_t i=0;i<3;i++) {
    msg_Out()<<"    p["<<i<<"]    = "<<split.Momentum(i)<<"\n";     
  }
    msg_Out()<<"    p[spec] = "<<split.SpectatorMomentum()<<"\n";     
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & tcut) :
  p_splitter(splitter), p_spectator(spectator), p_kernel(NULL), m_mode(splitting_mode::diff),
  m_Q2(0.), m_Q2red(m_Q2), m_y(0), m_eta(0),  
  m_tstart(t), m_tcut(tcut),
  m_isend(false), 
  p_weight(NULL),
  m_kinscheme(-1), m_clustered(0)
{
  for (size_t i=0;i<3;i++) {
    if (i<2)  m_t[i] = m_z[i] = m_phi[i] = 0.;
    if (i==0) m_t[i] = t;
    m_ztilde[i] = m_m2s[i] = 0.;
    m_outs[i]   = NULL;
  }
  if (p_splitter!=NULL && p_spectator!=NULL) {
    m_Q2   = ATOOLS::dabs(((p_splitter->Beam()>0?1.:-1.)*p_splitter->Mom() +
			   (p_spectator->Beam()>0?1.:-1.)*p_spectator->Mom()).Abs2());
    m_eta  = splitter->XB();
  }
  s_cnt++;
}

Splitting::~Splitting() {
  if (p_weight) { delete p_weight; p_weight = NULL; }
  s_cnt--;
}

void Splitting::InitSplitting(const ATOOLS::Mass_Selector * ms) {
  m_eta       = p_splitter->XB(); 
  m_msplit2   = sqr(ms->Mass(p_splitter->Flav()));
  m_mspect2   = sqr(ms->Mass(p_spectator->Flav()));
  m_Q2red     = m_Q2 - m_mspect2;
  m_ismassive = (m_mspect2>0.);
  if (p_kernel==0 || p_kernel->GetSF()==0) {
    msg_Error()<<METHOD<<" throws error: "<<p_kernel<<(p_kernel==0?0:p_kernel->GetSF())<<"\n";
    return;
  }
  m_nflavs = p_kernel->GetSF()->GetFlavs().size();
  for (size_t i=0;i<3;i++) {
    if (i<m_nflavs) {
      m_m2s[i]   = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[i]));
      m_Q2red   -= m_m2s[i];
      if (m_m2s[i]>0.) m_ismassive = true;
    }
    m_outs[i]  = NULL;
  }
}

bool Splitting::InitLimits() {
  double disc = 1.-4.*m_tcut/m_Q2;
  if (disc>=0.) {
    double delta = sqrt(disc);
    m_zmin = (1.-delta)/2.;
    m_zmax = (1.+delta)/2.;
    return true;
  }
  return false;
}



