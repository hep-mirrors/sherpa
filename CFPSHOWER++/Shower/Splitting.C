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

ostream & CFPSHOWER::operator<<(ostream &s,const Splitting & split) {
  s<<"Splitting(splitter = "<<split.GetSplitter()->Id()<<", "
   <<"spectator = "<<split.GetSpectator()->Id()<<", "
   <<"t = "<<split.T()<<", z = "<<split.Z()<<") for "
   <<(*split.GetKernel());
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & t0) :
  p_splitter(splitter), p_spectator(spectator), p_kernel(NULL), p_weight(NULL),
  m_t(t), m_t0(t0), m_Q2(-1.), m_z(-1.), m_phi(-1.), m_sijk(-1.),
  m_x(-1.), m_y(-1.), m_eta(splitter->XB()),
  m_mij2(0.), m_mi2(0.), m_mj2(0.), m_mk2(0.),
  m_specmom(p_spectator->Mom()), 
  m_kinscheme(-1), m_clustered(0)
{
  for (size_t i=0;i<3;i++) {
    m_moms[i] = ATOOLS::Vec4D(0.,0.,0.,0.);
    m_cols[i] = Color(0,0);
    p_outs[i] = NULL;
  }
  SetQ2();
  s_cnt++;
}

bool Splitting::InitKinematics(const ATOOLS::Mass_Selector * ms) {
  SetMasses(ms);
  m_Q2 = m_sijk = (p_splitter->Mom()+p_spectator->Mom()).Abs2();
  if (m_mi2>100. || m_mj2>100. ||
      p_kernel->GetSF()->GetFlavs()[1]==Flavour(kf_t) ||
      p_kernel->GetSF()->GetFlavs()[1]==Flavour(kf_t).Bar()) {}
  m_Q2  = dabs(m_sijk - m_mi2 - m_mj2 - m_mk2);
  m_eta = p_splitter->XB(); 
  //if (m_sijk<sqr(sqrt(m_mi2)+sqrt(m_mj2)+sqrt(m_mk2))) return false;
  return true;
}

void Splitting::SetMasses(const Mass_Selector * ms) {
  m_mij2 = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[0]));
  m_mi2  = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[1]));
  m_mj2  = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[2]));
  m_mk2  = sqr(ms->Mass(p_spectator->Flav()));
}


