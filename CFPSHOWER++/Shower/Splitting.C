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
   <<"t = "<<split.t();
  for (size_t i=0;i<3;i++) s<<", z["<<i<<"] = "<<split.z(i)<<") for "
   <<(*split.GetKernel());
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & t0) :
  p_splitter(splitter), p_spectator(spectator), p_kernel(NULL), 
  m_t(t), m_tstart(t), m_t0(t0),
  m_Q2(-1.), m_s(-1.), m_x(-1.), m_y(-1.), m_eta(splitter->XB()), m_phi(-1.), 
  p_weight(NULL),
  m_kinscheme(-1), m_clustered(0)
{
  m_msplit2 = (p_splitter->Mom()).Abs2();
  m_specmom = p_spectator->Mom();
  m_mspect2 = m_specmom.Abs2();
  for (size_t i=0;i<3;i++) {
    m_m2[i]   = m_z[i] = m_sij[i] = 0.;
    m_moms[i] = ATOOLS::Vec4D(0.,0.,0.,0.);
    m_cols[i] = Color(0,0);
    p_outs[i] = NULL;
  }
  Set_Q2();
  s_cnt++;
}

void Splitting::InitKinematics(const ATOOLS::Mass_Selector * ms) {
  SetMasses(ms);
  m_s   = (p_splitter->Mom()+p_spectator->Mom()).Abs2();
  m_Q2  = dabs(m_s - m_m2[0] - m_m2[1] - m_m2[2] - m_mspect2);
  m_eta = p_splitter->XB(); 
}

void Splitting::UpdateSpectatorMomentum() {
  p_spectator->SetMom(m_specmom);
  // update Bjorken-x variables of particles, where appropriate
  if (p_kernel->GetType()==kernel_type::IF ||
      p_kernel->GetType()==kernel_type::II) p_splitter->SetXB();
  if (p_kernel->GetType()==kernel_type::FI ||
      p_kernel->GetType()==kernel_type::II) p_spectator->SetXB();
}


void Splitting::SetMasses(const Mass_Selector * ms) {
  if (p_kernel==0) return;
  m_msplit2 = sqr(ms->Mass(p_splitter->Flav()));
  m_mspect2 = sqr(ms->Mass(p_spectator->Flav()));
  for (size_t i=0;i<p_kernel->GetSF()->GetFlavs().size();i++)
    m_m2[i] = sqr(ms->Mass(p_kernel->GetSF()->GetFlavs()[i]));
}


