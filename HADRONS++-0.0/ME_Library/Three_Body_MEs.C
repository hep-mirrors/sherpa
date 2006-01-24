#include "Three_Body_MEs.H"
#include "Message.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

//##############################################################################
//##############################################################################


P_3P_Dalitz::P_3P_Dalitz(int _nout, ATOOLS::Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs), m_pseudo_1(-1), m_pseudo_2(-1), m_pseudo_3(-1),
  m_const(1.),m_liny(0.),m_linx(0.),m_linyphase(0.),m_linxphase(0.),
  m_quady(0.),m_quadx(0.),m_quadyphase(0.),m_quadxphase(0.),
  m_phaseliny(0.),m_phaselinx(0.),m_phasequady(0.),m_phasequadx(0.)
{
  m_metype = string("P_3P_Dalitz");
  for (int i=1;i<4;i++) {
    if( p_flavs[i].IsAnti() != p_flavs[0].IsAnti() ) {
      m_pseudo_3 = i;
      m_mProd3   = p_flavs[i].Mass();
    }
    else if( m_pseudo_1 == -1 ) {
      m_pseudo_1 = i;
      m_mProd1   = p_flavs[i].Mass();
    }
    else {
      m_pseudo_2 = i;
      m_mProd2   = p_flavs[i].Mass();
    }
  }
  m_mDecayer = p_flavs[0].Mass();
}

void P_3P_Dalitz::SetModelParameters( GeneralModel _md )
{
  m_const = _md("const",1.);
  m_liny  = _md("liny",0.); 
  m_linx = _md("linx",0.); 
  m_quady = _md("quady",0.); 
  m_quadx = _md("quadx",0.);
  
  m_linyphase  = _md("linyphase",0.); 
  m_linxphase  = _md("linxphase",0.);
  m_quadyphase = _md("quadyphase",0.); 
  m_quadxphase = _md("quadxphase",0.);
  m_phaseliny  = _md("phaseliny",0.); 
  m_phaselinx  = _md("phaselinx",0.);
  m_phasequady = _md("phaseliny",0.); 
  m_phasequadx = _md("phaselinx",0.);
}

double P_3P_Dalitz::operator()( const Vec4D * _p )
{
  // kinematic variables
  double s1 = (_p[0]-_p[m_pseudo_1]).Abs2();
  double s2 = (_p[0]-_p[m_pseudo_2]).Abs2();
  double s3 = (_p[0]-_p[m_pseudo_3]).Abs2();
  double s0 = (s1+s2+s3)/3.;
  double x  = (s1-s2)/s0;
  double y  = (s3-s0)/s0;

  // amplitude
  Complex ampl = Complex(m_const,0.);
  ampl += y*Complex(m_liny + m_linyphase*cos(m_phaseliny),m_linyphase*sin(m_phaseliny));
  ampl += x*Complex(m_linx + m_linxphase*cos(m_phaselinx),m_linyphase*sin(m_phaselinx));
  ampl += sqr(y)*Complex(m_quady + m_quadyphase*cos(m_phasequady),m_quadyphase*sin(m_phasequady));
  ampl += sqr(x)*Complex(m_quadx + m_quadxphase*cos(m_phasequadx),m_quadyphase*sin(m_phasequadx));

  // amplitude squared
  double ampl_sq = std::abs(ampl*conj(ampl));

  return ampl_sq;
}


//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################

P_GammaFF::P_GammaFF(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs), m_phot(-1), m_f1(-1), m_f2(-1)
{ 
  m_metype = string("P_GammaFF");
  for (int i=1;i<4;i++) {
    if (p_flavs[i].IsPhoton()) m_phot=i;
    else { 
      if (m_f1<0) m_f1 = i;
             else m_f2 = i;
    } 
  }
}

double P_GammaFF::operator()(const Vec4D * p)
{
  double pref = 4.*M_PI/137.;
  
  Vec4D  q2    = p[m_f1]+p[m_f2];
  double q22   = q2.Abs2();
  double pfpfb = p[m_f1]*p[m_f2];
  double mfmfb = p_masses[m_f1]*p_masses[m_f2];
  double q1q2  = p[m_phot]*q2;
  double q1pf  = p[m_phot]*p[m_f1];
  double q1pfb = p[m_phot]*p[m_f2];
  
  return 4.*pref*((2.*pfpfb+4.*mfmfb)*sqr(q1q2) - 2.*q22*q1pf*q1pfb)/sqr(q22);
}

//##############################################################################
//##############################################################################
//##############################################################################

P_2PGamma::P_2PGamma(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs), m_phot(-1), m_p1(-1), m_p2(-1)
{ 
  m_metype = string("P_2PGamma");
  for (int i=1;i<4;i++) {
    if (p_flavs[i].IsPhoton()) m_phot=i;
    else { 
      if (m_p1<0) m_p1 = i;
             else m_p2 = i;
    } 
  }
}

void   P_2PGamma::operator()( 
    const ATOOLS::Vec4D  * _p, 
    std::vector<Complex> * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  _ampls_tensor->clear();
  Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
  _ampls_tensor->push_back( ampl );
  _indices->clear();
}

double P_2PGamma::operator()(const Vec4D * p)
{
  double kp1  = p[m_phot]*p[m_p1];
  double kp2  = p[m_phot]*p[m_p2];
  double p1p2 = p[m_p1]*p[m_p2];
  return 
    -(p_masses2[m_p1]*sqr(kp2)+p_masses2[m_p2]*sqr(kp1)-2.*kp1*kp2*p1p2)/pow(p_masses[0],6.);
}

//##############################################################################
//##############################################################################
//##############################################################################

P_P2Gamma::P_P2Gamma(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs), m_phot1(-1), m_phot2(-1), m_p(-1)
{ 
  m_metype = string("P_P2Gamma");
  for (int i=1;i<4;i++) {
    if (!p_flavs[i].IsPhoton()) m_p=i;
    else { 
      if (m_phot1<0) m_phot1 = i;
                else m_phot2 = i;
    } 
  }
  m_mrho2 = sqr(Flavour(kf::rho_770).Mass());
  m_grho2 = sqr(Flavour(kf::rho_770).Width());

//  cout<<"New P2Gamma "<<p_flavs[0]<<" -> "
//      <<p_flavs[1]<<" "<<p_flavs[2]<<" "<<p_flavs[3]
//      <<" "<<m_phot1<<m_phot2<<m_p<<endl; 
}

void   P_P2Gamma::operator()( 
    const ATOOLS::Vec4D  * _p, 
    std::vector<Complex> * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  _ampls_tensor->clear();
  Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
  _ampls_tensor->push_back( ampl );
  _indices->clear();
}

double P_P2Gamma::operator()(const Vec4D * p)
{
  double s = (p[m_phot1]+p[m_phot2]).Abs2();
  double t = (p[m_p]+p[m_phot2]).Abs2();
  double u = (p[m_p]+p[m_phot1]).Abs2();
 
  Complex aV = 
    (t+p_masses2[0])/(sqr(t-m_mrho2)+m_mrho2*m_grho2)*
    Complex(t-m_mrho2,-sqrt(m_mrho2*m_grho2)) +
    (u+p_masses2[0])/(sqr(u-m_mrho2)+m_mrho2*m_grho2)*
    Complex(u-m_mrho2,-sqrt(m_mrho2*m_grho2));
  Complex bV = 
    1./(sqr(t-m_mrho2)+m_mrho2*m_grho2)*
    Complex(t-m_mrho2,-sqrt(m_mrho2*m_grho2)) +
    1./(sqr(u-m_mrho2)+m_mrho2*m_grho2)*
    Complex(u-m_mrho2,-sqrt(m_mrho2*m_grho2));
  
  double kkp   = p[m_phot1]*p[m_phot2];
  double pk    = p[0]*p[m_phot1];
  double pkp   = p[0]*p[m_phot2];
  return abs(
    (sqr(kkp)) * real(aV*conj(aV)) +
    (-4.*kkp*pk*pkp-sqr(kkp)*p_masses2[0]) * 2.* real(aV*conj(bV)) +
    (2.*sqr(pk)*sqr(pkp)-2.*kkp*p_masses2[0]*pk*pkp+sqr(kkp)*sqr(p_masses2[0])) * real(bV*conj(bV)) );
}

//##############################################################################
//##############################################################################
//##############################################################################

P_3P_DalitzDef::P_3P_DalitzDef(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs),
  m_allpions(true), m_allsame(true),
  m_pi0(1),m_pim(2),m_pip(3)
{   
  m_metype = string("P_3P_DalitzDef");
  for (int i=1;i<4;i++) {
    if (p_flavs[i]!=Flavour(kf::pi)) {
      if (p_flavs[i]!=Flavour(kf::pi_plus)) {
	m_allpions=false;
	break;
      }
      else m_allsame = false;
    }
  }
  if (!m_allsame) {
    for (int i=1;i<4;i++) {
      if (p_flavs[i]==Flavour(kf::pi))                              m_pi0 = i;
      if (p_flavs[i]==Flavour(kf::pi_plus) && !p_flavs[i].IsAnti()) m_pip = i;
      if (p_flavs[i]==Flavour(kf::pi_plus) &&  p_flavs[i].IsAnti()) m_pim = i;
    }
  }
//  cout<<"New DalitzDef "<<p_flavs[0]<<" -> "
//      <<p_flavs[1]<<" "<<p_flavs[2]<<" "<<p_flavs[3]
//      <<" "<<m_allpions<<m_allsame<<" "<<m_pi0<<m_pip<<m_pim<<endl; 
}

double P_3P_DalitzDef::operator()(const Vec4D * p)
{
  double s0  = p_masses[0]/3.+p_masses[1]+p_masses[2]+p_masses[3];
  double mpi = (p_masses[1]+p_masses[2]+p_masses[3])/3.;
  if (m_allpions) {
    double A = 0.;
    A += 1.+3.*((p[m_pip]+p[m_pim]).Abs2()-s0)/(p_masses[0]-mpi);
    if (!m_allsame) return sqr(A);
    A += 1.+3.*((p[m_pip]+p[m_pi0]).Abs2()-s0)/(p_masses[0]-mpi);
    A += 1.+3.*((p[m_pi0]+p[m_pim]).Abs2()-s0)/(p_masses[0]-mpi);    
    return sqr(A);
  }
  return 1.;
}

//##############################################################################
//##############################################################################
//##############################################################################

//P_3P_Dalitz::P_3P_Dalitz(int _nout,Flavour * _flavs) :
//  HD_ME_Base(_nout,_flavs),
//  m_a1(0.), m_a2(0.), m_b1(0.), m_b2(0.), m_c(0.)
//{ 
//  m_metype = string("P_3P_Dalitz");
//  cout<<"New Dalitz "<<p_flavs[0]<<" -> "
//      <<p_flavs[1]<<" "<<p_flavs[2]<<" "<<p_flavs[3]<<endl;
//}
//
//double P_3P_Dalitz::operator()(const Vec4D * moms)
//{
//  double s_a = (moms[0]-moms[1]).Abs2();
//  double s_b = (moms[0]-moms[2]).Abs2();
//  double s_c = (moms[0]-moms[3]).Abs2();
//  double Q   = p_masses[0]-p_masses[1]-p_masses[2]-p_masses[3];
//  double x   = sqrt(3.)/(2.*p_masses[0]*Q)*(s_c-s_b);
//  double y   = 3./(2.*p_masses[0]*Q)*(sqr(p_masses[0]-p_masses[1])-s_a)-1.;
//  
//  return (1. + m_a1*y + m_a2*y*y + m_b1*x + m_b2*x*x + m_c*x*y);
//}
//
//void P_3P_Dalitz::SetDalitzParameters(std::vector<double> & _dals)
//{
//  if (_dals.size()!=5) {
//    msg.Error()<<"Error in P_3P_Dalitz::SetDalitzParameters : "<<endl
//	       <<"   Not enough parameters, only "<<_dals.size()<<","<<endl
//	       <<"   Will ignore this."<<endl;
//    m_a1=m_a2=m_b1=m_b2=m_c = 0.;
//  }
//}
//
