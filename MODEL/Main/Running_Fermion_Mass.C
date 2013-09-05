#include "MODEL/Main/Running_Fermion_Mass.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;

const double ZETA3 = 1.2020569031595942855;
const double ZETA5 = 1.036927755143369926331;

Running_Fermion_Mass::Running_Fermion_Mass(ATOOLS::Flavour _flav,double _yukmass,
					   Running_AlphaS * _as) :
  p_as(as), m_fl(_flav)
{
  if (m_fl.Mass(true)<1.0||(!_flav.IsQuark())||_yukmass<1.0) {
    m_order = 0;
    p_as    = NULL;
    return;
  }
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  m_runbelowpole = dataread.GetValue<int>("RUN_MASS_BELOW_POLE",0);
  if (m_runbelowpole)
    msg_Debugging()<<METHOD<<"(): "<<m_fl<<" mass runs below pole."<<std::endl;
  m_polemass = GetMSBarMass(m_fl.Mass(true));
  msg_Tracking()<<METHOD<<":("<<m_fl<<") m_{pole} = "<<m_fl.Mass(true)
		<<" -> m_{MSbar} = "<<m_polemass<<".\n";
  m_type    = std::string("Running Mass");
  m_name    = "Mass_"+ToString(m_fl);
  m_defval  = m_polemass;
  m_a       = (*p_as)(sqr(m_polemass));
  m_order = p_as->Order()+1;
}

double Running_Fermion_Mass::GetMSBarMass(const double &mp) const
{
  double a=(*p_as)(mp*mp)/M_PI;
  double C=1.0-4.0/3.0*a;
  if (p_as->Order()>0) {
    double nl=m_fl.Kfcode()-1.0, nh=1.0;
    double d12=7.0/128.0-3.0/4.0*ZETA3+1.0/2.0*sqr(M_PI)*log(2.0)-5.0/16.0*sqr(M_PI);
    double d22=-1111.0/384.0+3.0/8.0*ZETA3-1.0/4.0*sqr(M_PI)*log(2.0)+1.0/12.0*sqr(M_PI);
    double d32=71.0/96.0+1.0/12.0*sqr(M_PI);
    double d42=143.0/96.0-1.0/6.0*sqr(M_PI);
    C+=4.0/3.0*a*a*(4.0/3.0*d12+3.0*d22+0.5*nl*d32+0.5*nh*d42);
    if (p_as->Order()>1) {
      double a4=0.5174790616738994;
      double d13=-2969.0/768.0-1.0/16.0*sqr(M_PI)*ZETA3-81.0/16.0*ZETA3+5.0/8.0*ZETA5+29.0/4.0*sqr(M_PI)*log(2.0)
	+1.0/2.0*sqr(M_PI*log(2.0))-613.0/192.0*sqr(M_PI)-1.0/48.0*pow(M_PI,4.0)-1.0/2.0*pow(log(2.0),4.0)-12.0*a4;
      double d23=13189.0/4608.0-19.0/16.0*sqr(M_PI)*ZETA3-773.0/96.0*ZETA3+45.0/16.0*ZETA5-31.0/72.0*sqr(M_PI)*log(2.0)
	-31.0/36.0*sqr(M_PI*log(2.0))+509.0/576.0*sqr(M_PI)+65.0/432.0*pow(M_PI,4.0)-1.0/18.0*pow(log(2.0),4.0)-4.0/3.0*a4;
      double d33=-1322545.0/124416.0+51.0/64.0*sqr(M_PI)*ZETA3+1343.0/288.0*ZETA3-65.0/32.0*ZETA5
	-115.0/72.0*sqr(M_PI)*log(2.0)+11.0/36.0*sqr(M_PI*log(2.0))-1955.0/3456*sqr(M_PI)-179.0/3456.0*pow(M_PI,4.0)
	+11.0/72.0*pow(log(2.0),4.0)+11.0/3.0*a4;
      double d43=1283.0/576.0+55.0/24.0*ZETA3-11.0/9.0*sqr(M_PI)*log(2.0)+2.0/9.0*sqr(M_PI*log(2.0))+13.0/18.0*sqr(M_PI)
	-119.0/2160.0*pow(M_PI,4.0)+1.0/9.0*pow(log(2.0),4.0)+8.0/3.0*a4;
      double d53=1067.0/576.0-53.0/24.0*ZETA3+8.0/9.0*sqr(M_PI)*log(2.0)-1.0/9.0*sqr(M_PI*log(2.0))-85.0/108.0*sqr(M_PI)
	+91.0/2160.0*pow(M_PI,4.0)+1.0/9.0*pow(log(2.0),4.0)+8.0/3.0*a4;
      double d63=70763.0/15552.0+89.0/144.0*ZETA3+11.0/18.0*sqr(M_PI)*log(2.0)-1.0/9.0*sqr(M_PI*log(2.0))
	+175.0/432.0*sqr(M_PI)+19.0/2160.0*pow(M_PI,4.0)-1.0/18.0*pow(log(2.0),4.0)-4.0/3.0*a4;
      double d73=144959.0/15552.0+1.0/8.0*sqr(M_PI)*ZETA3-109.0/144.0*ZETA3-5.0/8.0*ZETA5+32.0/9.0*sqr(M_PI)*log(2.0)
	+1.0/18.0*sqr(M_PI*log(2.0))-449.0/144.0*sqr(M_PI)-43.0/1080.0*pow(M_PI,4.0)-1.0/18.0*pow(log(2.0),4.0)-4.0/3.0*a4;
      double d83=-5917.0/3888.0+2.0/9.0*ZETA3+13.0/108.0*sqr(M_PI);
      double d93=-9481.0/7776.0+11.0/18.0*ZETA3+4.0/135.0*sqr(M_PI);
      double d103=-2353.0/7776.0-7.0/18.0*ZETA3-13.0/108.0*sqr(M_PI);
      C+=4.0/3.0*a*a*a*
	(sqr(4.0/3.0)*d13+4.0/3.0*3.0*d23+sqr(3.0)*d33
	 +4.0/3.0*0.5*nl*d43+4.0/3.0*0.5*nh*d53+3.0*0.5*nl*d63
	 +3.0*0.5*nh*d73+sqr(0.5)*nl*nh*d83+sqr(0.5)*nh*d93
	 +sqr(0.5*nl)*d103);
    }
  }
  return mp*C;
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
