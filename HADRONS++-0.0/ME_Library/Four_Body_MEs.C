#include "Four_Body_MEs.H"
#include "Message.H"
#include "Tools.H"
#include "XYZFuncs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Higgs_TauTau_2Pi2Neutrino::Higgs_TauTau_2Pi2Neutrino( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl)
{
  m_metype = string("Higgs_TauTau_2Pi2Neutrino");
}

void Higgs_TauTau_2Pi2Neutrino::SetModelParameters( GeneralModel _md )
{
  m_pseudo = int(0.5 + _md("pseudo",0.0)); // scalar or pseudoscalar
}

void Higgs_TauTau_2Pi2Neutrino::operator()(
                       const Vec4D         * _p,
                       Amplitude_Tensor    * amps,
                       int                   k0_n )
{
  // Decay: B+ --> pi+ pi- nu_tau nu_taub
  //        0  --> 1   2   3      4
  Complex i(0.0,1.0);

  // XYZ functions for "higgs -> tau tau" part
  Vec4D   p_higgs[3];
  Flavour f_higgs[3];
  p_higgs[0]  = _p[0];
  f_higgs[0]  = p_flavs[0];
  p_higgs[1]  = _p[2]+_p[3];
  f_higgs[1]  = Flavour(kf::tau);
  p_higgs[2]  = _p[1]+_p[4];
  f_higgs[2]  = Flavour(kf::tau).Bar();
  XYZFunc F_higgs(2,&p_higgs[0],&f_higgs[0],k0_n);

  // XYZ functions for "tau- -> pi- nu_tau" part
  Vec4D   p_tau[3];
  Flavour f_tau[3];
  p_tau[0]  = _p[2]+_p[3];
  f_tau[0]  = Flavour(kf::tau);
  p_tau[1]  = _p[2];
  f_tau[1]  = p_flavs[2];
  p_tau[2]  = _p[3];
  f_tau[2]  = p_flavs[3];
  XYZFunc F_tau(2,&p_tau[0],&f_tau[0],k0_n);

  // XYZ functions for "tau+ -> pi+ nu_taubar" part
  Vec4D   p_taubar[3];
  Flavour f_taubar[3];
  p_taubar[0]  = _p[1]+_p[4];
  f_taubar[0]  = Flavour(kf::tau).Bar();
  p_taubar[1]  = _p[1];
  f_taubar[1]  = p_flavs[1];
  p_taubar[2]  = _p[4];
  f_taubar[2]  = p_flavs[4];
  XYZFunc F_taubar(2,&p_taubar[0],&f_taubar[0],k0_n);

  // propagator denominators
  double tauq2 =    (_p[2]+_p[3]).Abs2();
  double taubarq2 = (_p[1]+_p[4]).Abs2();
  double mtau =     Flavour(kf::tau).Mass();
  double Gtau =     Flavour(kf::tau).Width();
  Complex prop1 = i/(tauq2    - sqr(mtau) + i*mtau*Gtau);
  Complex prop2 = i/(taubarq2 - sqr(mtau) + i*mtau*Gtau);
  Complex props = prop1*prop2;

  Complex X1, X2, Y;
  for( int hnutau1=0; hnutau1<2; hnutau1++ ) {
    for( int hnutau2=0; hnutau2<2; hnutau2++ ) {
      Complex ampl(0.0,0.0);
      for( int htau1=0; htau1<2; htau1++ ) {
        for( int htau2=0; htau2<2; htau2++ ) {
          X1 = F_tau.X(2,hnutau1,1,0,htau1,0.0,1.0);
          X2 = F_taubar.X(0,htau2,1,2,hnutau2,0.0,1.0);

          if(m_pseudo) Y  = F_higgs.Y(1,htau1,2,htau2,1.0,-1.0);
          else         Y  = -i*F_higgs.Y(1,htau1,2,htau2,1.0,1.0);
          
          ampl += Y*X1*X2;
        }
      }
      ampl = ampl*props;

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,0));
      spins.push_back(make_pair(1,0));
      spins.push_back(make_pair(2,0));
      spins.push_back(make_pair(3,hnutau1));
      spins.push_back(make_pair(4,hnutau2));
      amps->InsertAmplitude(ampl,spins);
    }
  }
  F_higgs.Delete();
  F_tau.Delete();
  F_taubar.Delete();
}
