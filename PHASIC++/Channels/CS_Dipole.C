#include "PHASIC++/Channels/CS_Dipole.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Flow.H"

using namespace ATOOLS;
using namespace PHASIC;

CS_Dipole::CS_Dipole(const ATOOLS::Cluster_Leg &lij,
		     const ATOOLS::Cluster_Leg &lk,
		     const ATOOLS::Cluster_Leg &li,
		     const ATOOLS::Cluster_Leg &lj):
  m_lij(lij), m_lk(lk), m_li(li), m_lj(lj),
  p_vegas(NULL), m_alpha(1.0), m_oldalpha(1.0),
  m_amin(0.0), m_amax(1.0)
{
  Reset();
}

CS_Dipole::~CS_Dipole() 
{
  if (p_vegas) delete p_vegas;
}

void CS_Dipole::InitVegas(Process_Base *const proc)
{
  m_id=ToString(ID(m_lij.Id()).front())+"_"
    +ToString(ID(m_lk.Id()).front())+"__"
    +ToString(ID(m_li.Id()).front())+"_"
    +ToString(ID(m_lj.Id()).front());
  p_vegas = new Vegas(3,100,proc->Name()+"_"+m_id,0);
  // p_vegas->SetCheckMode(0);
}

void CS_Dipole::AddPoint(const double &weight)
{
  if (m_weight==0.0) {
    p_vegas->AddPoint(0.0,m_rn);
    return;
  }
  double wgt(m_weight!=0.0?sqr(weight)/m_weight:0.0);
  ++m_np;
  m_sum+=wgt;
  m_sum2+=sqr(wgt);
  m_max=ATOOLS::Max(m_max,dabs(wgt));
  p_vegas->AddPoint(weight,m_rn);
}

void CS_Dipole::Map(Cluster_Amplitude *const ampl)
{
  size_t i(ID(m_li.Id()).front());
  size_t j(ID(m_lj.Id()).front());
  size_t ij(ID(m_lij.Id()).front());
  if (i>ij) 
    for (size_t n(ij);n<i;++n) ampl->SwapLegs(n,n+1);
  else if (ij>i)
    for (size_t n(ij);n>i;--n) ampl->SwapLegs(n,n-1);
  for (size_t n(ampl->Legs().size()-1);n>j;--n)
    ampl->SwapLegs(n,n-1);
  msg_Debugging()<<*ampl<<"\n";
  for (size_t n(0);n<ampl->Legs().size();++n)
    ampl->Leg(n)->SetId(1<<n);
}

void CS_Dipole::ReMap(Cluster_Amplitude *const ampl)
{
  size_t i(ID(m_li.Id()).front());
  size_t j(ID(m_lj.Id()).front());
  size_t ij(ID(m_lij.Id()).front());
  for (size_t n(j);n<ampl->Legs().size()-1;++n)
    ampl->SwapLegs(n,n+1);
  if (i>ij) 
    for (size_t n(i);n>ij;--n) ampl->SwapLegs(n,n-1);
  else if (ij>i)
    for (size_t n(i);n<ij;++n) ampl->SwapLegs(n,n+1);
  for (size_t n(0);n<ampl->Legs().size();++n)
    ampl->Leg(n)->SetId(1<<n);
}

void CS_Dipole::Reset()
{
  m_np=0.0;
  m_sum=m_sum2=m_max=0.0;
}

double CS_Dipole::Lambda
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

double CS_Dipole::Phi
(Vec4D pijt,Vec4D pkt,Vec4D pi,const bool ii) const
{
  Vec4D ktt(0.0,cross(Vec3D(pijt),Vec3D(pkt)));
  Poincare cms(pijt+pkt);
  cms.Boost(pijt);
  cms.Boost(pi);
  Poincare zax(pijt,Vec4D::ZVEC);
  if (!ii && ktt.PSpat2()>1.0e-6) zax.Rotate(ktt);
  else ktt=Vec4D(0.0,1.0,1.0,0.0);
  zax.Rotate(pi);
  Poincare xax(ktt,Vec4D::XVEC);
  xax.Rotate(pi);
  return pi.Phi();
}

std::ostream &PHASIC::operator<<(std::ostream &ostr,const CS_Dipole &dip) 
{
  return ostr<<"("<<&dip<<"): m_a = "<<dip.Alpha()<<" {\n"
	     <<"  "<<dip.LIJ()<<"\n  "<<dip.LK()<<"\n}";
}

