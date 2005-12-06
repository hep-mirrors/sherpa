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
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == kf::e ||
        p_flavs[i].Kfcode() == kf::mu )
      { m_lep = i; break; }     // that's the lepton
  }
  // find the corresponding neutrino
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == p_flavs[m_lep].Kfcode()+1 ) m_nulep = i;
  }
}

void K_Meson_Lepton::SetModelParameters( GeneralModel _md )
{
  m_GF2 = sqr(_md("GF",1.16639e-5));
  m_cL = (0.,_md("F_K", 1.0));
  m_cR = (0.0,0.0);
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
  m_nulep(-1)
{
  m_metype = string("K_Meson_Lepton");
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == kf::e ||
        p_flavs[i].Kfcode() == kf::mu )
    { m_lep = i; break; }     // that's the lepton
  }
  // find the corresponding neutrino
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == p_flavs[m_lep].Kfcode()+1 ) m_nulep = i;
  }
}

void K_Meson_SemiLeptonic::SetModelParameters( GeneralModel _md )
{
  m_GF2 = sqr(_md("GF",1.16639e-5));
  m_cL = (0.,_md("F_K", 1.0));
  m_cR = (0.0,0.0);
}

double K_Meson_SemiLeptonic::Using_Hels( const Vec4D * _p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {
    ret += norm( F.X(m_nulep, 0, m_lep, h, m_cR, m_cL) );
  }
  F.Delete();
  return ret*0.5;
}

double K_Meson_SemiLeptonic::operator()( const Vec4D *_p )
{
  double T (1.);
  T = Using_Hels(_p);
  return T*m_GF2;
}
