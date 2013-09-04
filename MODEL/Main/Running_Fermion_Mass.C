#include "MODEL/Main/Running_Fermion_Mass.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;

const double ZETA3 = 1.2020569031595942855;

Running_Fermion_Mass::Running_Fermion_Mass(ATOOLS::Flavour _flav,double _polemass,
					   Running_AlphaS * _as) :
  m_polemass(_polemass), p_as(as)
{
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  m_runbelowpole = dataread.GetValue<int>("RUN_MASS_BELOW_POLE",0);
  if (m_runbelowpole)
    msg_Debugging()<<METHOD<<"(): "<<_flav<<" mass runs below pole."<<std::endl;
  m_type    = std::string("Running Mass");
  m_name    = "Mass_"+ToString(_flav);
  m_defval  = m_polemass;
  m_a       = (*p_as)(sqr(m_polemass));
  if (_flav.Mass(true)<1.||(!_flav.IsQuark())||m_polemass<1.) {
    m_order = 0;
    p_as    = NULL;
    return;
  }
  m_order = p_as->Order()+1;
}

double Running_Fermion_Mass::Beta0(const double &nf) const
{
  return 1.0/4.0*(11.0-2.0/3.0*nf);
}

double Running_Fermion_Mass::Beta1(const double &nf) const
{
  return 1.0/16.0*(102.0-38.0/3.0*nf);
}

double Running_Fermion_Mass::Beta2(const double &nf) const
{
  return 1.0/64.0*(2857.0/2.0-5033.0/18.0*nf+325.0/54.0*nf*nf);
}

double Running_Fermion_Mass::Gamma0(const double &nf) const
{
  return 1.0;
}

double Running_Fermion_Mass::Gamma1(const double &nf) const
{
  return 1/16.0*(202.0/3.0-20.0/9.0*nf);
}

double Running_Fermion_Mass::Gamma2(const double &nf) const
{
  return 1/64.0*(1249.0+(-2216.0/27.0-160.0/3.0*ZETA3)*nf-140.0/81.0*nf*nf);
}

double Running_Fermion_Mass::Series(const double &a,const int nf) const
{
  double s=1.0, c0=Gamma0(nf), b0=Beta0(nf);
  if (p_as->Order()>0) {
    double b1=Beta1(nf), c1=Gamma1(nf);
    double A1=-b1*c0/sqr(b0)+c1/b0;
    s+=a*A1;
    if (p_as->Order()>1) {
      double b2=Beta2(nf), c2=Gamma2(nf);
      double A2=c0/sqr(b0)*(sqr(b1)/b0-b2)-b1*c1/sqr(b0)+c2/b0;
      s+=a*a/2.0*(A1*A1+A2);
    }
  }
  return pow(a/b0/2.,c0/b0)*s;  
}

double Running_Fermion_Mass::operator()(double t) {
  if (m_order==0) return m_polemass;
  if (t<0.) t = -t;
  if (!m_runbelowpole && t<sqr(m_polemass)) return m_polemass;
  double nf=p_as->Nf(t);
  return m_polemass/Series(m_a,nf)*Series((*p_as)(t),nf);
}

void Running_Fermion_Mass::SelfTest() {
  double m_test = m_polemass/2.;
  for (int i=0;i<100;i++) {
    m_test += m_polemass/20.*i;
    std::cout<<"  "<<m_test<<" "<<(*this)(sqr(m_test))<<std::endl;
  }
}
