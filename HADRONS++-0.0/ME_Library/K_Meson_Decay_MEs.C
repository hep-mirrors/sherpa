#include "K_Meson_Decay_MEs.H"
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

K_Meson_Lepton::K_Meson_Lepton( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_lep(-1),
  m_nulep(-1)
{
  m_metype = string("K_Meson_Lepton");
  if( p_flavs[1].Kfcode() == kf::e || p_flavs[1].Kfcode() == kf::mu ) { 
    m_lep=1;
    m_nulep=2; 
  }
  else { m_lep=2; m_nulep=1; }
}

void K_Meson_Lepton::SetModelParameters( GeneralModel _md )
{
  m_GF2 = sqr(_md("GF",1.16639e-5));
  m_cL = (0.0,_md("F_K", 1.0));  m_cR = (0.0,0.0);
}

double K_Meson_Lepton::Using_Hels( const Vec4D * _p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {
      ret += norm( F.X(m_nulep, 0, m_lep, h, m_cR, m_cL) );
  }
  F.Delete();
  return ret*0.5;
}

double K_Meson_Lepton::operator()( const Vec4D *_p )
{
  double T (1.);
  T = Using_Hels(_p);
  return T*m_GF2;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  semileptonic decay  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_Meson_SemiLeptonic::K_Meson_SemiLeptonic( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_lep(-1),
  m_nulep(-1),
  m_pi(-1)
{
  m_metype = string("K_Meson_SemiLeptonic");
  // find the lepton-, neutrino- and pion-indices
  for( int i=1; i<4; i++ ) {
    switch (p_flavs[i].Kfcode()) {
      case kf::e:   case kf::mu:    m_lep = i;   break;
      case kf::nue: case kf::numu:  m_nulep = i; break;
      case kf::pi:                  m_pi = i;    break;
    }
  }
}

void K_Meson_SemiLeptonic::SetModelParameters( GeneralModel _md )
{
  m_GF2 = sqr(_md("GF",1.16639e-5));
  double f_plus = _md("f_plus",1.0); //only temporarily.
  double f_minus = _md("f_minus",1.0); //only temporarily.
  m_cL_K = (0.0,(f_plus+f_minus));  m_cR_K = (0.0,0.0);
  m_cL_pi = (0.0,(f_plus-f_minus));  m_cR_pi = (0.0,0.0);
}

double K_Meson_SemiLeptonic::Using_Hels( const Vec4D * _p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<8; h++ ) {
    ret += norm( F.X(m_nulep, 0, m_lep, h, m_cR_K, m_cL_K) 
               + F.X(m_nulep, m_pi, m_lep, h, m_cR_pi, m_cL_pi) );
  }
  F.Delete();
  return ret*0.5;
}

double K_Meson_SemiLeptonic::operator()( const Vec4D *_p )
{
  double T (1.);
  T = Using_Hels(_p);
  return T*m_GF2*0.25;
}

