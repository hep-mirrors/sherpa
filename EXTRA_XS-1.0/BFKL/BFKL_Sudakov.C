#include "BFKL_Sudakov.H"

#include "LO_LEV.H"
#include "Random.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace ATOOLS;

size_t BFKL_Sudakov::s_nf(5);

BFKL_Sudakov::BFKL_Sudakov():
  p_levs(new LEV_Group()),
  m_weight(0.0), m_ktexp(-1.0) {}

BFKL_Sudakov::~BFKL_Sudakov()
{
  delete p_levs;
}

bool BFKL_Sudakov::Initialize()
{
  if (m_splitmode&1) p_levs->Add(new G_GG());
  for (size_t i(1);i<=s_nf;++i) {
    if (m_splitmode&(1<<i)) {
      Flavour fl((kf::code)i);
      p_levs->Add(new Q_GQ(fl));
      p_levs->Add(new Q_GQ(fl.Bar()));
      p_levs->Add(new Q_QG(fl));
      p_levs->Add(new Q_QG(fl.Bar()));
      p_levs->Add(new G_QQ(fl));
      p_levs->Add(new G_QQ(fl.Bar()));
    }
  }
  return true;
}

double BFKL_Sudakov::DicePolynomial(const double &xmin,const double &xmax,
				    const double &a,const double &rn) const
{
  if (a==-1.0) return xmin*pow(xmax/xmin,rn);
  return pow(pow(xmax,a+1.0)*rn+(1.0-rn)*pow(xmin,a+1.0),1.0/(a+1.0));
}

double BFKL_Sudakov::WeightPolynomial(const double &xmin,const double &xmax,
				      const double &a,const double &x) const
{
  if (a==-1.0) return log(xmax/xmin)*x;
  return (pow(xmax,a+1.0)-pow(xmin,a+1.0))/(a+1.0)*pow(x,-a);
}

bool BFKL_Sudakov::Init()
{
  m_asmax=(*MODEL::as)(m_kt2min);
  p_levs->SetKT2Min(m_kt2min);
  m_y=m_ya;
  return true;
}

bool BFKL_Sudakov::Dice()
{
  if (m_kt2max<=m_kt2min) return false;
  double rn[4];
  for (short unsigned int i(0);i<4;++i) rn[i]=ran.Get();
  // kt limit
  p_levs->SetKT2Max(m_kt2max);
  // get last integral
  m_integral=p_levs->MajorIntegral(m_flq);
  // set branching probability
  double qt2(m_kt2min);
  switch (m_ktscheme) {
  case 0:
    break;
  case 1: {
    qt2=m_q.PPerp2();
    break;
  }
  default:
    THROW(fatal_error,"Invalid scale scheme.");
  }
  double b0(MODEL::as->Beta0(qt2)/M_PI);
  double l2(qt2*exp(-1.0/(b0*(*MODEL::as)(qt2))));
  m_gamma=1.0/(2.0*M_PI*b0)*
    log(log(m_q.PPerp2()/l2)/log(m_kt2min/l2))*m_integral;
  // dice new y
  if (m_ya>m_yb) m_y+=log(rn[0])/m_gamma;
  else m_y-=log(rn[0])/m_gamma;
  // dice new kt2
  m_kt2=l2*exp(DicePolynomial(log(m_kt2min/l2),
			      log(m_kt2max/l2),m_ktexp,rn[1]));
  // dice new phi
  m_phi=2.0*M_PI*rn[2];
  // select splitting
  p_levs->SelectLEV(rn[3]);
  return m_ya>m_yb?m_y>m_yb:m_y<m_yb;
}

bool BFKL_Sudakov::Approve(const ATOOLS::Vec4D &k1,
			   const ATOOLS::Vec4D &q1,
			   const ATOOLS::Vec4D &k2,
			   const ATOOLS::Vec4D &q2)
{
  // lower qt and kt limits
  if (q2.PPerp2()<m_kt2min || 
      k2.PPerp2()<m_kt2min) return false;
  // rejection weight
  double weight(1.0);
  // splitting weight
  weight*=(m_splitweight=p_levs->Selected()->Value(k1,q1,k2,q2)/
	   p_levs->Selected()->MajorValue(k1,q1,k2,q2));
  // veto emission
  if (weight<ran.Get()) return false;
  // reweighting weight
  m_weight=1.0;
  // divide by \frac{d}{dy}\ln\Delta(y,\bar y)
  m_weight/=m_gamma;
  // set scale
  double kt2(m_kt2min);
  switch (m_ktscheme) {
  case 0:
    break;
  case 1: {
    kt2=k2.PPerp2();
    break;
  }
  default:
    THROW(fatal_error,"Invalid scale scheme.");
  }
  // coupling & splitting weight
  m_weight*=(*MODEL::as)(kt2)/(2.0*M_PI*m_kt2)*m_integral;
  double qt2(q1.PPerp2()), b0(MODEL::as->Beta0(qt2)/M_PI);
  double l2(qt2*exp(-1.0/(b0*(*MODEL::as)(qt2))));
  // kt integration domain weight
  m_weight*=m_kt2*WeightPolynomial
    (log(m_kt2min/l2),log(m_kt2max/l2),m_ktexp,log(m_kt2/l2));
  // mass of t-channel particle
  double mq12(sqr(p_levs->Selected()->GetA().PSMass()));
  // propagator weight for massive quarks
  m_weight*=kt2/(kt2+mq12);
  return true;
}

