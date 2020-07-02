#include "CFPSHOWER++/Calculators/SF_Base.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::SF_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

std::ostream & CFPSHOWER::operator<<(std::ostream &s,const subtract::code & sub) {
  if      (sub==subtract::none) s<<"none";
  else if (sub==subtract::coll) s<<"coll";
  else if (sub==subtract::soft) s<<"soft";
  else if (sub==subtract::both) s<<"both";
  return s;
}

SF_Base::SF_Base(const Kernel_Info & info) :
  m_split(info.GetSplit()), m_nout(info.GetFlavs().size()),
  m_flavs(info.GetFlavs()), m_tags(info.TagSequence()),
  m_ismassive(false),
  m_logtype(info.LogType()), m_type(info.Type()),
  m_name("generic SF"),
  m_CMW(info.KFactor()), m_softcorr(info.SoftCorrection()), m_endpoint(info.Endpoint()),
  m_subtract(subtract::none),
  m_weight(0.)
{
  m_m.resize(m_nout);
  m_m2.resize(m_nout);
  m_moms.resize(m_nout+1);
  m_pp.resize(m_nout+1);
  for (size_t i=0;i<m_nout+1;i++) m_pp[i].resize(m_nout+1);
}  

bool SF_Base:: Init(const Splitting & split,const ATOOLS::Mass_Selector * msel) {
  m_spect  = split.GetSpectator()->Flav();
  m_msplit = msel->Mass(m_split); m_msplit2 = sqr(m_msplit);
  m_mspect = msel->Mass(m_spect); m_mspect2 = sqr(m_mspect);
  if (m_mspect>0.) m_ismassive = true;
  for (size_t i=0;i<m_nout;i++) {
    m_m[i]  = msel->Mass(m_flavs[i]);
    if (m_m[i]<1.e-6) m_m[i] = m_m2[i] = 0.; 
    else m_m2[i] = sqr(m_m[i]);
    if (m_m[i]>0.) m_ismassive = true;
  }
  m_psplit = split.GetSplitter()->Mom();
  m_pspect = split.GetSpectator()->Mom();
  m_pboth  = m_psplit + m_pspect;
  m_shat   = m_pboth.Abs2();  m_Ehat = sqrt(m_shat);
  return true;
}

Vec4D SF_Base::SumMomenta(const Configuration & config) {
  Vec4D momsum(0.,0.,0.,0.);
  for (Parton_List::const_iterator pit=config.begin();pit!=config.end();pit++) {
    momsum += (*pit)->Mom();
  }
  return momsum;
}

Vec4D SF_Base::LT(const Vec4D &a,const Vec4D &b,const Vec4D &c) {
  double b1c2 = b[1]*c[2]-b[2]*c[1], b0c2 = b[0]*c[2]-b[2]*c[0];
  double b2c3 = b[2]*c[3]-b[3]*c[2], b0c1 = b[0]*c[1]-b[1]*c[0];
  double b3c1 = b[3]*c[1]-b[1]*c[3], b3c0 = b[3]*c[0]-b[0]*c[3];
  
  double t =  a[1]*b2c3 + a[2]*b3c1 + a[3]*b1c2;
  double x = -a[0]*b2c3 - a[2]*b3c0 - a[3]*b0c2;
  double y = -a[0]*b3c1 + a[1]*b3c0 + a[3]*b0c1;
  double z = -a[0]*b1c2 + a[1]*b0c2 - a[2]*b0c1;
  return Vec4D(t,-x,-y,-z);
}

double SF_Base::Lambda(const double & a,const double & b,const double & c) {
  double lambda2 = Lambda2(a,b,c); 
  if (lambda2<0.) {
    msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<b<<", "<<c<<") yields nan.\n"
	       <<"   return 0. and hope for the best.\n";
    return 0.;
  }
  return sqrt(lambda2);
}

double SF_Base::Lambda2(const double & a,const double & b,const double & c) {
  return sqr(a-b-c)-4.*b*c;
}
