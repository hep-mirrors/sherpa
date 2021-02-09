#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CFPSHOWER {
  size_t Splitting::s_cnt=0;
}

ostream & CFPSHOWER::operator<<(ostream &s,Splitting & split) {
  s<<"Splitting(splitter = "<<split.GetSplitter()->Id()<<", "
   <<"spectator = "<<split.GetSpectator()->Id()<<", "
   <<"t = "<<split.T()<<", ";
  if (split.GetKernel()->KinScheme()==kin_type::code::CS)
    s<<"z = "<<split.Z()<<", ";
  else if (split.GetKernel()->KinScheme()==kin_type::code::PanGlobal)
    s<<"y = "<<split.Y()<<", ";
  s<<"    kernel = "<<(*split.GetKernel());
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & tcut) :
  p_splitter(splitter), p_spectator(spectator),
  p_kernel(NULL), p_weight(NULL),
  m_Q2(t), m_tstart(t), m_tcut(tcut), m_t(t),
  m_y(0.), m_eta(0.), m_phi(0.), m_zmin(0.), m_zmax(1.), m_z(0.), m_kt2(tcut),
  m_momscale(tcut)
{
  if (splitter!=NULL && spectator!=NULL) {
    m_Q2   = dabs(((splitter->Beam()>0?1.:-1.)  * splitter->Mom() +
		   (spectator->Beam()>0?1.:-1.) * spectator->Mom()).Abs2());
  }
  s_cnt++;
}

Splitting::~Splitting() {
  if (p_weight) { delete p_weight; p_weight = NULL; }
  s_cnt--;
}

bool Splitting::InitLimits() {
  m_moms.reserve(p_kernel->GetFlavs().size()+1);
  if (p_kernel->KinScheme()==kin_type::PanGlobal) {
    return true;
  }
  else if (p_kernel->KinScheme()==kin_type::CS) {
    double disc = 1.-4.*m_tcut/m_Q2;
    if (disc>=0.) {
      double delta = sqrt(disc);
      m_zmin = (1.-delta)/2.;
      m_zmax = (1.+delta)/2.;
      return true;
    }
  }
  return false;
}








