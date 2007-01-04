#include "Two_Body_MEs.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

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
  m_global  = GF*SQRT_05*(p_flavs[m_boson].PSMass());
  m_cR  = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL  = Complex(0.,_md("a",1.)+_md("b",1.));
}
 
void F_VF::operator()( 
    const Vec4D         * _p,
    Spin_Amplitudes    * amps,
    int                   k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  Vec4D eps[2];  
  Complex ampl (0.,0.);
  for( int hboson=0; hboson<3; hboson++ ) {      // sum over lambda (s?,l,t1,t2)
    Tools::ComplexBosonPolarizationVector(_p[m_boson],hboson,eps);
    for(int hf1=0;hf1<2;hf1++) {
      for(int hf2=0;hf2<2;hf2++) {
        ampl = Complex(0.,0.);
        if( !eps[0].IsZero() && !eps[0].Nan() )
          ampl += F.X(m_fermion,hf2, eps[0], 0,hf1, m_cR, m_cL);
        if( !eps[1].IsZero() && !eps[1].Nan() )
          ampl += Complex(0.,1.)*F.X(m_fermion,hf2, eps[1], 0,hf1, m_cR, m_cL);
        ampl *= m_global;
        
        vector<pair<int,int> > spins;
        spins.push_back(make_pair(0,hf1));
        spins.push_back(make_pair(m_fermion,hf2));
        spins.push_back(make_pair(m_boson,hboson));
        amps->Insert(ampl,spins);
      }
    }
  }
  F.Delete();
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
}

void V_FF::operator()(
    const Vec4D         * _p,
    Spin_Amplitudes    * amps,
    int                   k0_n)
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  Vec4D eps[2];
  Complex ampl(0.,0.);
  for( int hboson=0; hboson<3; hboson++ ) {      // sum over lambda (s?,l,t1,t2)
    Tools::ComplexBosonPolarizationVector(_p[0],hboson,eps);
    for(int hf1=0;hf1<2;hf1++) {
      for(int hf2=0;hf2<2;hf2++) {
        ampl = Complex(0.,0.);
        if( !eps[0].IsZero() && !eps[0].Nan() )
          ampl += F.X(m_ferm1,hf1, eps[0], m_ferm2,hf2, m_cR, m_cL);
        if( !eps[1].IsZero() && !eps[1].Nan() )
          ampl += Complex(0.,1.)*F.X(m_ferm1,hf1, eps[1], m_ferm2,hf2, m_cR, m_cL);
        ampl *= m_global;
        
        vector<pair<int,int> > spins;
        spins.push_back(make_pair(0,hboson));
        spins.push_back(make_pair(m_ferm1,hf1));
        spins.push_back(make_pair(m_ferm2,hf2));
        amps->Insert(ampl,spins);
      }
    }
  }
  F.Delete();
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_PP::V_PP( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl)
{
  m_metype = string("Vector_PseudoscalarPseudoscalar");
}

void V_PP::SetModelParameters( GeneralModel _md )
{
}

void V_PP::operator()(
                       const Vec4D         * _p,
                       Spin_Amplitudes    * amps,
                       int                   k0_n )
{
  Vec4D q = _p[1]-_p[2];
  Vec4D eps[2];
  Complex ampl (0.,0.);
  for( int hvector=0; hvector<3; hvector++ ) {      // sum over lambda (s,l,t1,t2)
    Tools::ComplexBosonPolarizationVector(_p[0],hvector,eps);

    vector<pair<int,int> > spins;
    spins.push_back(make_pair(0,hvector));
    spins.push_back(make_pair(1,0));
    spins.push_back(make_pair(2,0));

    ampl = Complex(0.,0.);
    ampl += eps[0]*q;
    ampl += Complex(0.0,1.0)*(eps[1]*q) ;

    amps->Insert(ampl,spins);
  }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_FF::S_FF( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl)
{
  m_metype = string("S_FF");
}

void S_FF::SetModelParameters( GeneralModel _md )
{
  m_pseudo = int(0.5 + _md("pseudo",0.0)); // scalar or pseudoscalar
}

void S_FF::operator()(
                       const Vec4D         * _p,
                       Spin_Amplitudes    * amps,
                       int                   k0_n )
{
  Complex ampl(0.,0.);
  Complex i(0.0,1.0);
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  for( int htau1=0; htau1<2; htau1++ ) {
    for( int htau2=0; htau2<2; htau2++ ) {
      if(m_pseudo) ampl = F.Y(1,htau1,2,htau2,1.0,-1.0);
      else         ampl = -i*F.Y(1,htau1,2,htau2,1.0,1.0);
      
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,0));
      spins.push_back(make_pair(1,htau1));
      spins.push_back(make_pair(2,htau2));
      amps->Insert(ampl,spins);
    }
  }
  F.Delete();
}
