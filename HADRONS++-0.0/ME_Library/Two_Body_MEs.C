#include "Two_Body_MEs.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

P_2Gamma::P_2Gamma(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs)
{ 
  m_metype = string("P_2Gamma");
}

void   P_2Gamma::operator()( 
    const ATOOLS::Vec4D  * _p, 
    std::vector<Complex> * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  _ampls_tensor->clear();
  _ampls_tensor->push_back( p_masses2[0]/2. );
  _indices->clear();
}

double P_2Gamma::operator()(const Vec4D * moms)
{
  return sqr(p_masses2[0])/4.;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  top -> bottom + W         %%%%%%%%%%%%%%%%%%%%%%%%%%%
//  tau -> rho + neutrino     %%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_VF::F_VF( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_fermion(-1),
  m_boson(-1)
{
  m_metype = string("Fermion_VectorFermion");
  for( int i=1; i<3; i++ ) {
    if( p_flavs[i].IsFermion() )  m_fermion = i;        // that's the fermion
    else                          m_boson  = i;         // that's the boson
  }
}
 
void F_VF::SetModelParameters( GeneralModel _md ) 
{ 
  double GF = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
  m_global  = GF*SQRT_05*Flavour(kf::W).PSMass();
  m_cR  = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL  = Complex(0.,_md("a",1.)+_md("b",1.));

  m_MB  = Flavour(kf::W).PSMass();
  m_GB  = Flavour(kf::W).Width(); 
  m_MB2 = sqr(m_MB);
}
 
void F_VF::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  Vec4D  p  = _p[m_boson];
  double p2 = p.Abs2();
  _ampls_tensor->clear();
  Complex ampl (0.,0.);
  for( int h=0; h<4; h++ ) {        // for all hel. comb. (b,t)
    for( int l=0; l<4; l++ ) {      // sum over lambda (s,l,t1,t2)
      ATOOLS::Vec4D* eps = Tools::ComplexBosonPolarizationVector(p,l,sqr(80.419));
      ampl = Complex(0.,0.);
      if( !eps[0].IsZero() && !eps[0].Nan() ) ampl += F.X(m_fermion, eps[0], 0, h, m_cR, m_cL);
      if( !eps[1].IsZero() && !eps[1].Nan() ) ampl += Complex(0.,1.)*F.X(m_fermion, eps[1], 0, h, m_cR, m_cL);
      delete[] eps;
      ampl *= m_global;
      _ampls_tensor->push_back( ampl );
    }
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2 3)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(m_boson,3) );      // it has 4 polarisations !
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_fermion,1) );
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  W   -> lepton + neutrino  %%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
V_FF::V_FF( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_ferm1(-1),
  m_ferm2(-1)
{
  m_metype = string("Vector_FermionFermion");
  for( int i=1; i<3; i++ ) {
    if( p_flavs[i].IsAnti()==p_flavs[0].IsAnti() )  m_ferm1 = i;       
    else                                            m_ferm2 = i;
  }
}
 
void V_FF::SetModelParameters( GeneralModel _md ) 
{ 
  double GF = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
  m_global  = GF*SQRT_05*Flavour(kf::W).PSMass();
  m_cR  = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL  = Complex(0.,_md("a",1.)+_md("b",1.));

  m_MB  = Flavour(kf::W).PSMass();
  m_GB  = Flavour(kf::W).Width(); 
  m_MB2 = sqr(m_MB);
}
 
void V_FF::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  Vec4D  p  = _p[0];
  double p2 = p.Abs2();
  _ampls_tensor->clear();
  bool zero_eps;
  Complex ampl(0.,0.);
  for( int h=0; h<4; h++ ) {        // for all hel. comb. (b,t)
    for( int l=0; l<4; l++ ) {      // sum over lambda (s,l,t1,t2)
      ATOOLS::Vec4D* eps = Tools::ComplexBosonPolarizationVector(p,l,sqr(80.419));
      ampl = Complex(0.,0.);
      if( !eps[0].IsZero() && !eps[0].Nan() ) ampl += F.X(m_ferm1, eps[0], m_ferm2, h, m_cR, m_cL);
      if( !eps[1].IsZero() && !eps[1].Nan() ) ampl += Complex(0.,1.)*F.X(m_ferm1, eps[1], m_ferm2, h, m_cR, m_cL);
      delete[] eps;
      ampl *= m_global;
      _ampls_tensor->push_back( ampl );
    }
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2 3)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,3) );      // it has 4 polarisations !
  _indices->push_back( pair<int,int>(m_ferm2,1) );
  _indices->push_back( pair<int,int>(m_ferm1,1) );
}

