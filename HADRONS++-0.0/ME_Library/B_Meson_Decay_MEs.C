#include "B_Meson_Decay_MEs.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Histogram.H"
#include "Traces.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  leptonic decay  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Semileptonic_B_Meson::Semileptonic_B_Meson( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_had(-1),m_nu(-1),m_lep(-1)
{
  m_metype = string("Semileptonic_B_Meson");
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == kf::e ||
	p_flavs[i].Kfcode() == kf::mu ||
	p_flavs[i].Kfcode() == kf::tau ) m_lep = i; // that's the lepton
    else if( p_flavs[i].Kfcode() == kf::nue ||
	     p_flavs[i].Kfcode() == kf::numu ||
	     p_flavs[i].Kfcode() == kf::nutau ) m_nu = i; // that's the lepton
    else m_had = i;
  }
  m_Bmass   = p_flavs[0].Mass();
  m_hadmass = p_flavs[m_had].Mass();
}
 
void Semileptonic_B_Meson::SetModelParameters( GeneralModel _md ) 
{ 
}
 
double Semileptonic_B_Meson::operator()( const Vec4D *_p )
{
  double  T  = Using_Hels(_p);
  return T;
}

double Semileptonic_B_Meson::Using_Hels( const Vec4D * _p )
{
  Complex FFplus  = Formfactor_plus(_p);
  Complex FFminus = Formfactor_minus(_p);
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {				// only nulep and lep have a helicity !
    ret += norm(( F.X(m_lep,0,m_nu,h,0.,1.)/m_Bmass +
		  F.X(m_lep,m_had,m_nu,h,0.,1.)/m_hadmass)*FFplus +
		( F.X(m_lep,0,m_nu,h,0.,1.)/m_Bmass - 
		  F.X(m_lep,m_had,m_nu,h,0.,1.)/m_hadmass)*FFminus);
  }
  F.Delete();
  return ret;
} 
 
Complex Semileptonic_B_Meson::Formfactor_plus( const Vec4D *_p )
{
  return Complex(0.,0.);
}

Complex Semileptonic_B_Meson::Formfactor_minus( const Vec4D *_p )
{
  double vvprime = _p[0]*_p[m_had]/(m_Bmass*m_hadmass);
  double xi      = 1.-0.7*(vvprime-1.);
  return Complex(xi,0.);
}
