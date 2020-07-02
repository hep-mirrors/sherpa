#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
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

ostream & CFPSHOWER::operator<<(ostream &s,Splitting & split) {
  s<<"Splitting(splitter = "<<split.GetSplitter()->Id()<<", "
   <<"spectator = "<<split.GetSpectator()->Id()<<", "
   <<"t1 = "<<split.T(0)<<", "
   <<"t2 = "<<split.T(1)<<", "
   <<"z1 = "<<split.Z(0)<<", "
   <<"z2 = "<<split.Z(1)<<", "
   <<"    kernel = "<<(*split.GetKernel());
  //for (size_t i=0;i<3;i++) {
  //  msg_Out()<<"    p["<<i<<"]    = "<<split.GetKinematics()->m_moms[i+1]<<"\n";     
  //}
  //msg_Out()<<"    p[spec] = "<<split.GetKinematics()->m_moms[0]<<"\n";     
  return s;
}

Splitting::Splitting(Parton * splitter,Parton * spectator,
		     const double  & t, const double  & tcut) :
  p_splitter(splitter), p_spectator(spectator), p_kernel(NULL), 
  m_tstart(t), m_tcut(tcut), m_zmin(0.), m_zmax(1.),
  m_isendpoint(false), m_isclustered(false), 
  m_kinscheme(-1), 
  p_weight(NULL)
{
  for (size_t i=0;i<3;i++) {
    if (i<2) m_t[i] = m_z[i] = m_phi[i] = 0.;
    p_outs[i] = NULL;
  }
  m_t[0] = m_tstart;
  if (splitter!=NULL && spectator!=NULL) {
    m_Q2   = dabs(((splitter->Beam()>0?1.:-1.)  * splitter->Mom() +
		   (spectator->Beam()>0?1.:-1.) * spectator->Mom()).Abs2());
    m_xB = splitter->XB();
  }
  s_cnt++;
}

Splitting::~Splitting() {
  if (p_weight) { delete p_weight; p_weight = NULL; }
  s_cnt--;
}

bool Splitting::InitLimits() {
  double disc = 1.-4.*m_tcut/m_Q2;
  if (disc>=0.) {
    double delta = sqrt(disc);
    m_zmin = (1.-delta)/2.;
    m_zmax = (1.+delta)/2.;
    return true;
  }
  //msg_Out()<<METHOD<<" yields false: "<<m_tcut<<" / "<<m_Q2<<"\n";
  return false;
}








