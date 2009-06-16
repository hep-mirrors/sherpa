#include "HADRONS++/ME_Library/Three_Body_MEs.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "ATOOLS/Org/Message.H"
#include "HELICITIES/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

void B_Bpi_pwave::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = -m_cR;
}

void B_Bpi_pwave::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  vector<pair<int,int> > spins(3);
  spins[2] = make_pair(p_i[2],0);
  for(int h0=0; h0<2; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      Complex amp=F.X(1,h1, p[p_i[2]], 0,h0, m_cR, m_cL);
      amps->Insert(amp,spins);
    }
  }
}

DEFINE_ME_GETTER(B_Bpi_pwave,B_Bpi_pwave_Getter,"B_Bpi_pwave")

void B_Bpi_pwave_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B3_Bpi_pwave::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = m_cR;
}

void B3_Bpi_pwave::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  vector<pair<int,int> > spins(3);
  spins[2] = make_pair(p_i[2],0);
  for(int h0=0; h0<4; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      amps->Insert(p[p_i[2]]*F.Y13(1,h1, 0,h0, m_cR, m_cL),spins);
    }
  }
}

DEFINE_ME_GETTER(B3_Bpi_pwave,B3_Bpi_pwave_Getter,"B3_Bpi_pwave")

void B3_Bpi_pwave_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B_Bpi_swave::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = m_cR;
}

void B_Bpi_swave::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  vector<pair<int,int> > spins(3);
  spins[2] = make_pair(p_i[2],0);
  for(int h0=0; h0<2; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      amps->Insert(F.Y(1,h1, 0,h0, m_cR, m_cL),spins);
    }
  }
}

DEFINE_ME_GETTER(B_Bpi_swave,B_Bpi_swave_Getter,"B_Bpi_swave")

void B_Bpi_swave_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B3_Bpi_dwave::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = -m_cR;
}

void B3_Bpi_dwave::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D q=p[p_i[2]];
  vector<pair<int,int> > spins(3);
  spins[2] = make_pair(p_i[2],0);
  for(int h0=0; h0<4; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      Complex amp=q*F.X13(1,h1,q,0,h0,m_cR,m_cL);
      amps->Insert(-amp, spins);
    }
  }
}

DEFINE_ME_GETTER(B3_Bpi_dwave,B3_Bpi_dwave_Getter,"B3_Bpi_dwave")

void B3_Bpi_dwave_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B_Bphoton_M1::SetModelParameters( GeneralModel _md )
{
  m_fac = -_md("f",1.0);
}

void B_Bphoton_M1::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D q=p[p_i[2]];
  vector<pair<int,int> > spins(3);

  Polarization_Vector eps(q,m_flavs[p_i[2]].HadMass());
  for(int h0=0; h0<2; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      for(int h2=0; h2<2; ++h2) {
        spins[2] = make_pair(p_i[2],h2);
        amps->Insert(m_fac*F.G(1,h1,conj(eps[h2]),0,h0),spins);
      }
    }
  }
}

DEFINE_ME_GETTER(B_Bphoton_M1,B_Bphoton_M1_Getter,"B_Bphoton_M1")

void B_Bphoton_M1_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B3_Bphoton_M1::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = m_cR;
}

void B3_Bphoton_M1::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D q=p[p_i[2]];
  vector<pair<int,int> > spins(3);

  Polarization_Vector eps(q,m_flavs[p_i[2]].HadMass());
  for(int h0=0; h0<4; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      for(int h2=0; h2<2; ++h2) {
        spins[2] = make_pair(p_i[2],h2);
        Vec4C cr=cross(conj(eps[h2]),p[p_i[0]],q);
        amps->Insert(cr*F.Y31(1,h1,0,h0,m_cR,m_cL),spins);
      }
    }
  }
}

DEFINE_ME_GETTER(B3_Bphoton_M1,B3_Bphoton_M1_Getter,"B3_Bphoton_M1")

void B3_Bphoton_M1_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B3_Bphoton_M1_2::SetModelParameters( GeneralModel _md )
{
  m_cR = Complex(0.0,_md("f",1.0));
  m_cL = -m_cR;
}

void B3_Bphoton_M1_2::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D q=p[p_i[2]];
  vector<pair<int,int> > spins(3);

  Polarization_Vector eps(q,m_flavs[p_i[2]].HadMass());
  for(int h0=0; h0<4; ++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for(int h1=0; h1<2; ++h1) {
      spins[1] = make_pair(p_i[1],h1);
      for(int h2=0; h2<2; ++h2) {
        spins[2] = make_pair(p_i[2],h2);
        Complex amp=conj(eps[h2])*F.X31(1,h1,q,0,h0,m_cR,m_cL);
        amp-=q*F.X31(1,h1,conj(eps[h2]),0,h0,m_cR,m_cL);
        amps->Insert(amp,spins);
      }
    }
  }
}

DEFINE_ME_GETTER(B3_Bphoton_M1_2,B3_Bphoton_M1_2_Getter,"B3_Bphoton_M1_2")

void B3_Bphoton_M1_2_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B_Bphoton_E1::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = -m_cR;
}

void B_Bphoton_E1::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D q=p[p_i[2]];
  vector<pair<int,int> > spins(3);

  double p0q=p[p_i[0]]*q;

  Polarization_Vector eps(q,m_flavs[p_i[2]].HadMass());
  for(int h2=0; h2<2; ++h2) {
    spins[2] = make_pair(p_i[2],h2);
    Vec4C epsstar=conj(eps[h2]);
    Complex p0eps=p[p_i[0]]*epsstar;
    for(int h0=0; h0<2; ++h0) {
      spins[0] = make_pair(p_i[0],h0);
      for(int h1=0; h1<2; ++h1) {
        spins[1] = make_pair(p_i[1],h1);
        Complex amp=p0q*F.X(1,h1,epsstar,0,h0,m_cR,m_cL)-
          p0eps*F.X(1,h1,q,0,h0,m_cR,m_cL);
        amps->Insert(amp,spins);
      }
    }
  }
}

DEFINE_ME_GETTER(B_Bphoton_E1,B_Bphoton_E1_Getter,"B_Bphoton_E1")

void B_Bphoton_E1_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void B3_Bphoton_E1::SetModelParameters( GeneralModel _md )
{
  m_cR = _md("f",1.0);
  m_cL = m_cR;
}

void B3_Bphoton_E1::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  XYZFunc F(m_n,p,m_flavs,Tools::k0,m_anti,p_i);

  Vec4D p2=p[p_i[2]];
  vector<pair<int,int> > spins(3);

  double p0p2=p[p_i[0]]*p2;

  Polarization_Vector eps(p2,m_flavs[p_i[2]].HadMass());
  for(int h2=0; h2<2; ++h2) {
    spins[2] = make_pair(p_i[2],h2);
    Vec4C epsstar=conj(eps[h2]);
    Complex p0eps=p[p_i[0]]*epsstar;
    for(int h0=0; h0<4; ++h0) {
      spins[0] = make_pair(p_i[0],h0);
      for(int h1=0; h1<2; ++h1) {
        spins[1] = make_pair(p_i[1],h1);
        Complex amp=p0p2*(epsstar*F.Y31(1,h1,0,h0,m_cR,m_cL))-
          p0eps*(p2*F.Y31(1,h1,0,h0,m_cR,m_cL));
        amps->Insert(amp,spins);
      }
    }
  }
}

DEFINE_ME_GETTER(B3_Bphoton_E1,B3_Bphoton_E1_Getter,"B3_Bphoton_E1")

void B3_Bphoton_E1_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

void QQ_PGG::SetModelParameters( GeneralModel _md )
{
  m_min_mass2 = _md("min_mass2",1.0);
}

void QQ_PGG::operator()(
                      const Vec4D             * p,
                      HELICITIES::Spin_Amplitudes * amps)
{
  Vec4D mom_gg=p[p_i[2]]+p[p_i[3]];
  double mass2=mom_gg.Abs2();

  if (mass2<m_min_mass2)
    amps->CreateTrivial(Complex(0.0,0.0));
  else
    amps->CreateTrivial(Complex(1.0,0.0));
}

DEFINE_ME_GETTER(QQ_PGG,QQ_PGG_Getter,"QQ_PGG")

void QQ_PGG_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}

//##############################################################################
//##############################################################################

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

void   P_3P_Dalitz::operator()( 
    const ATOOLS::Vec4D  * _p, 
    HELICITIES::Spin_Amplitudes * amps)
{
  Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
  vector<pair<int,int> > spins;
  spins.push_back(make_pair(0,0));
  spins.push_back(make_pair(1,0));
  spins.push_back(make_pair(2,0));
  spins.push_back(make_pair(3,0));
  amps->Add(ampl,spins);
  
}

double P_3P_Dalitz::operator()( const Vec4D * _p )
{
  // kinematic variables
  double s1 = (_p[p_i[0]]-_p[p_i[1]]).Abs2();
  double s2 = (_p[p_i[0]]-_p[p_i[2]]).Abs2();
  double s3 = (_p[p_i[0]]-_p[p_i[3]]).Abs2();
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


// //##############################################################################
// //##############################################################################
// //##############################################################################
// //##############################################################################
// 
// P_GammaFF::P_GammaFF(int _nout,Flavour * _flavs) :
//   HD_ME_Base(_nout,_flavs), m_phot(-1), m_f1(-1), m_f2(-1)
// { 
//   m_metype = string("P_GammaFF");
//   for (int i=1;i<4;i++) {
//     if (p_flavs[i].IsPhoton()) m_phot=i;
//     else { 
//       if (m_f1<0) m_f1 = i;
//              else m_f2 = i;
//     } 
//   }
// }
//  
// void   P_GammaFF::operator()( 
//     const ATOOLS::Vec4D  * _p, 
//     std::vector<Complex> * _ampls_tensor, 
//     std::vector<std::pair<int,int> > * _indices,
//     int                    k0_n )
// {
//   _ampls_tensor->clear();
//   Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
//   _ampls_tensor->push_back( ampl );
//   _indices->clear();
// }
// 
// 
// double P_GammaFF::operator()(const Vec4D * p)
// {
//   double pref = 4.*M_PI/137.;
//   
//   Vec4D  q2    = p[m_f1]+p[m_f2];
//   double q22   = q2.Abs2();
//   double pfpfb = p[m_f1]*p[m_f2];
//   double mfmfb = p_masses[m_f1]*p_masses[m_f2];
//   double q1q2  = p[m_phot]*q2;
//   double q1pf  = p[m_phot]*p[m_f1];
//   double q1pfb = p[m_phot]*p[m_f2];
//   
//   return 4.*pref*((2.*pfpfb+4.*mfmfb)*sqr(q1q2) - 2.*q22*q1pf*q1pfb)/sqr(q22);
// }
// 
// //##############################################################################
// //##############################################################################
// //##############################################################################
// 
// P_2PGamma::P_2PGamma(int _nout,Flavour * _flavs) :
//   HD_ME_Base(_nout,_flavs), m_phot(-1), m_p1(-1), m_p2(-1)
// { 
//   m_metype = string("P_2PGamma");
//   for (int i=1;i<4;i++) {
//     if (p_flavs[i].IsPhoton()) m_phot=i;
//     else { 
//       if (m_p1<0) m_p1 = i;
//              else m_p2 = i;
//     } 
//   }
// }
// 
// void   P_2PGamma::operator()( 
//     const ATOOLS::Vec4D  * _p, 
//     std::vector<Complex> * _ampls_tensor, 
//     std::vector<std::pair<int,int> > * _indices,
//     int                    k0_n )
// {
//   _ampls_tensor->clear();
//   Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
//   _ampls_tensor->push_back( ampl );
//   _indices->clear();
// }
// 
// double P_2PGamma::operator()(const Vec4D * p)
// {
//   double kp1  = p[m_phot]*p[m_p1];
//   double kp2  = p[m_phot]*p[m_p2];
//   double p1p2 = p[m_p1]*p[m_p2];
//   return 
//     -(p_masses2[m_p1]*sqr(kp2)+p_masses2[m_p2]*sqr(kp1)-2.*kp1*kp2*p1p2)/pow(p_masses[0],6.);
// }
// 
// //##############################################################################
// //##############################################################################
// //##############################################################################
// 
// P_P2Gamma::P_P2Gamma(int _nout,Flavour * _flavs) :
//   HD_ME_Base(_nout,_flavs), m_phot1(-1), m_phot2(-1), m_p(-1)
// { 
//   m_metype = string("P_P2Gamma");
//   for (int i=1;i<4;i++) {
//     if (!p_flavs[i].IsPhoton()) m_p=i;
//     else { 
//       if (m_phot1<0) m_phot1 = i;
//                 else m_phot2 = i;
//     } 
//   }
//   m_mrho2 = sqr(Flavour(kf_rho_770).Mass());
//   m_grho2 = sqr(Flavour(kf_rho_770).Width());
// 
// //  cout<<"New P2Gamma "<<p_flavs[0]<<" -> "
// //      <<p_flavs[1]<<" "<<p_flavs[2]<<" "<<p_flavs[3]
// //      <<" "<<m_phot1<<m_phot2<<m_p<<endl; 
// }
// 
// void   P_P2Gamma::operator()( 
//     const ATOOLS::Vec4D  * _p, 
//     std::vector<Complex> * _ampls_tensor, 
//     std::vector<std::pair<int,int> > * _indices,
//     int                    k0_n )
// {
//   _ampls_tensor->clear();
//   Complex ampl = csqrt( (*this)(_p) );        // call uncorrelated
//   _ampls_tensor->push_back( ampl );
//   _indices->clear();
// }
// 
// double P_P2Gamma::operator()(const Vec4D * p)
// {
//   double s = (p[m_phot1]+p[m_phot2]).Abs2();
//   double t = (p[m_p]+p[m_phot2]).Abs2();
//   double u = (p[m_p]+p[m_phot1]).Abs2();
//  
//   Complex aV = 
//     (t+p_masses2[0])/(sqr(t-m_mrho2)+m_mrho2*m_grho2)*
//     Complex(t-m_mrho2,-sqrt(m_mrho2*m_grho2)) +
//     (u+p_masses2[0])/(sqr(u-m_mrho2)+m_mrho2*m_grho2)*
//     Complex(u-m_mrho2,-sqrt(m_mrho2*m_grho2));
//   Complex bV = 
//     1./(sqr(t-m_mrho2)+m_mrho2*m_grho2)*
//     Complex(t-m_mrho2,-sqrt(m_mrho2*m_grho2)) +
//     1./(sqr(u-m_mrho2)+m_mrho2*m_grho2)*
//     Complex(u-m_mrho2,-sqrt(m_mrho2*m_grho2));
//   
//   double kkp   = p[m_phot1]*p[m_phot2];
//   double pk    = p[0]*p[m_phot1];
//   double pkp   = p[0]*p[m_phot2];
//   return abs(
//     (sqr(kkp)) * real(aV*conj(aV)) +
//     (-4.*kkp*pk*pkp-sqr(kkp)*p_masses2[0]) * 2.* real(aV*conj(bV)) +
//     (2.*sqr(pk)*sqr(pkp)-2.*kkp*p_masses2[0]*pk*pkp+sqr(kkp)*sqr(p_masses2[0])) * real(bV*conj(bV)) );
// }
// 
// //##############################################################################
// //##############################################################################
// //##############################################################################
// 
// P_3P_DalitzDef::P_3P_DalitzDef(int _nout,Flavour * _flavs) :
//   HD_ME_Base(_nout,_flavs),
//   m_allpions(true), m_allsame(true),
//   m_pi0(1),m_pim(2),m_pip(3)
// {   
//   m_metype = string("P_3P_DalitzDef");
//   for (int i=1;i<4;i++) {
//     if (p_flavs[i]!=Flavour(kf_pi)) {
//       if (p_flavs[i]!=Flavour(kf_pi_plus)) {
// 	m_allpions=false;
// 	break;
//       }
//       else m_allsame = false;
//     }
//   }
//   if (!m_allsame) {
//     for (int i=1;i<4;i++) {
//       if (p_flavs[i]==Flavour(kf_pi))                              m_pi0 = i;
//       if (p_flavs[i]==Flavour(kf_pi_plus) && !p_flavs[i].IsAnti()) m_pip = i;
//       if (p_flavs[i]==Flavour(kf_pi_plus) &&  p_flavs[i].IsAnti()) m_pim = i;
//     }
//   }
// //  cout<<"New DalitzDef "<<p_flavs[0]<<" -> "
// //      <<p_flavs[1]<<" "<<p_flavs[2]<<" "<<p_flavs[3]
// //      <<" "<<m_allpions<<m_allsame<<" "<<m_pi0<<m_pip<<m_pim<<endl; 
// }
// 
// double P_3P_DalitzDef::operator()(const Vec4D * p)
// {
//   double s0  = p_masses[0]/3.+p_masses[1]+p_masses[2]+p_masses[3];
//   double mpi = (p_masses[1]+p_masses[2]+p_masses[3])/3.;
//   if (m_allpions) {
//     double A = 0.;
//     A += 1.+3.*((p[m_pip]+p[m_pim]).Abs2()-s0)/(p_masses[0]-mpi);
//     if (!m_allsame) return sqr(A);
//     A += 1.+3.*((p[m_pip]+p[m_pi0]).Abs2()-s0)/(p_masses[0]-mpi);
//     A += 1.+3.*((p[m_pi0]+p[m_pim]).Abs2()-s0)/(p_masses[0]-mpi);    
//     return sqr(A);
//   }
//   return 1.;
// }

DEFINE_ME_GETTER(P_3P_Dalitz,P_3P_Dalitz_Getter,"P_3P_Dalitz")

void P_3P_Dalitz_Getter::PrintInfo(std::ostream &str,const size_t width) const {
  str<<"No Information";
}
