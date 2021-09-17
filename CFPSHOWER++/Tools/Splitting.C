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
  if (split.GetKernel()->KinScheme()==kin_type::code::CataniSeymour)
    s<<"z = "<<split.Z()<<", ";
  else if (split.GetKernel()->KinScheme()==kin_type::code::Alaric)
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
  m_IntPhi(1.)
{
  if (splitter!=NULL && spectator!=NULL) {
    m_psplit  = p_splitter->Mom();
    m_pspect  = p_spectator->Mom();
    m_pdipole = m_psplit + m_pspect;
    m_Q2      = m_pdipole.Abs2();  m_Q = sqrt(m_Q2);
    //m_Q2      = dabs(((p_splitter->Beam()>0?1.:-1.)  * m_psplit +
    //		      (p_spectator->Beam()>0?1.:-1.) * m_pspect).Abs2());
  }
  s_cnt++;
}

Splitting::~Splitting() {
  if (p_weight) { delete p_weight; p_weight = NULL; }
  s_cnt--;
}

void Splitting::Init(Kernel * kernel,const ATOOLS::Mass_Selector * msel) {
  p_kernel    = kernel;
  m_ismassive = false;
  m_msplit = msel->Mass(p_splitter->Flav());  m_msplit2 = sqr(m_msplit);
  m_mspect = msel->Mass(p_spectator->Flav()); m_mspect2 = sqr(m_mspect);
  if (m_mspect>0.) m_ismassive = true;
  SF_Base * sf = p_kernel->GetSF();
  m_nout = sf->Flavs().size();
  m_moms.resize(m_nout+1);
  m_masses.resize(m_nout);
  m_masses2.resize(m_nout);
  m_masssum = 0.;
  for (size_t i=0;i<m_nout;i++) {
    m_masssum += m_masses[i]  = msel->Mass(sf->Flavs()[i]);
    if (m_masses[i]<1.e-6) m_masses[i] = m_masses2[i] = 0.; 
    else m_masses2[i] = sqr(m_masses[i]);
    if (m_masses[i]>0.) m_ismassive = true;
  }
}







