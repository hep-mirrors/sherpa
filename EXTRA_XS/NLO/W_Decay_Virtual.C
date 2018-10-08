#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Math/Poincare.H"
#define LOG_2 0.69314718055994530942
#include "EXTRA_XS/Main/ME2_Base.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define B_1(A,B,C,M)        PV_Bubble_1(A,B,C,M)
#define B_0p(A,B,C,M)       Master_Bubble_Prime(A,B,C,M)
#define B_1p(A,B,C,M)       PV_Bubble_1_Prime(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)
#define C_11(A,B,C,D,E,F,M) PV_Triangle_11(A,B,C,D,E,F,M)
#define C_12(A,B,C,D,E,F,M) PV_Triangle_12(A,B,C,D,E,F,M)
#define C_21(A,B,C,D,E,F,M) PV_Triangle_21(A,B,C,D,E,F,M)
#define C_22(A,B,C,D,E,F,M) PV_Triangle_22(A,B,C,D,E,F,M)
#define C_23(A,B,C,D,E,F,M) PV_Triangle_23(A,B,C,D,E,F,M)
#define C_24(A,B,C,D,E,F,M) PV_Triangle_24(A,B,C,D,E,F,M)



using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class W_Decay_Virtual : public ME2_Base {
    double m_fac;
    int                m_spins[4],m_ew, m_ew_scheme, m_qed, m_lo, m_own;
    Complex m_yqq, m_yll, m_cL_q, m_cL_l, m_cR_q, m_cR_l;
    Flavour    m_flavs[4];
    Flavour    m_quark_flavs[2];
    double m_masses[3];
    double m_m2, m_m2p, m_Qf, m_Qfp, m_If, m_Ifp,m_smod;
    double m_e, m_alpha, m_deltas, m_omega, m_mu2, chi, born;
    double chi1,chi2,qedterm,intterm,Zterm,m_kappa,qf,qe,vf,af,ve,ae,sin2tw,mass;
    double m_quark_masses[5];
    double MW2, GW2;
    double MZ2, GZ2;
    double MH2;
    double mt2;
    Complex muW2;
    Complex muZ2;
    Complex muH2;
    Complex mut2;
    Complex m_sw2_cms, m_sw_cms;
    Complex m_cw2_cms, m_cw_cms;
    Complex             m_sW, m_sW2, m_sW4;
    Complex             m_cW, m_cW2, m_cW4;
    ATOOLS::Vec4D m_moms[4];
    Complex    m_cL;
    Complex    m_cR;
    METOOLS::DivArrC One,Zero,m_res;
    ATOOLS::Vec4D m_p1,m_p2;
    double m_m_1,m_m_2,m_x1,m_x2,m_xx1,m_xx2;
    std::string widthscheme;
  public:
    W_Decay_Virtual(const Process_Info& pi, const Flavour_Vector& flavs);

    ~W_Decay_Virtual();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    Complex InfraredSubtractedME_0_0();
    METOOLS::DivArrC InfraredSubtractedME_0_1();
    Complex GetBeta_0_0();
    METOOLS::DivArrC GetBeta_0_1();
    METOOLS::DivArrC FAn(const double& s);
    METOOLS::DivArrC FZn(const double& s);
    METOOLS::DivArrC FZa(const double& s);
    METOOLS::DivArrC CT_L();
    METOOLS::DivArrC CT_R();
    METOOLS::DivArrC B_ff(const double& p2, const Complex& m12, const Complex& m22);
    // Self energies
    METOOLS::DivArrC Sigma_ZA(const double& p2);
    METOOLS::DivArrC Sigma_ZZ(const double& p2);
    METOOLS::DivArrC Sigma_W(const double& p2);
    METOOLS::DivArrC Sigma_H(const double& p2);
    // Derivatives of self energies
    METOOLS::DivArrC dSigma_ZZ(const double& p2);
    METOOLS::DivArrC dSigma_W(const double& p2);
    METOOLS::DivArrC dSigma_H(const double& p2);
    METOOLS::DivArrC dSigma_GamGam(const double& p2);
    METOOLS::DivArrC dSigma_GamGam0();

    // Fermion self energies
    METOOLS::DivArrC Sigma_ferm_R(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    METOOLS::DivArrC Sigma_ferm_L(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    METOOLS::DivArrC Sigma_ferm_S(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    // Derivatives of fermion self energies
    METOOLS::DivArrC dSigma_ferm_R(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);
    METOOLS::DivArrC dSigma_ferm_L(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);
    METOOLS::DivArrC dSigma_ferm_S(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);
    METOOLS::DivArrC dZW();
    METOOLS::DivArrC dZH();
    METOOLS::DivArrC dZZZ();
    METOOLS::DivArrC dZAA();
    METOOLS::DivArrC dZZA();
    METOOLS::DivArrC dZAZ();
    METOOLS::DivArrC dMW2();
    METOOLS::DivArrC dMZ2();
    METOOLS::DivArrC dMH2();
    METOOLS::DivArrC dZfermL(const double& m2, const double& m2p,
			     const double& Qf, const double& If);
    METOOLS::DivArrC dZfermR(const double& m2, const double& m2p,
			     const double& Qf, const double& If);
    METOOLS::DivArrC dm(const double& m2, const double& m2p,
			const double& Qf, const double& If);


    METOOLS::DivArrC dZe();
    METOOLS::DivArrC dcw();
    double YFS_Form_Factor();
    double IntP1();
    double IntE();
    double IntP2();
    double G(double);
    double IntG();
    double CalculateBeta(const ATOOLS::Vec4D&);
    void Print_Ren_Constants_Finite();
  };
}

using namespace METOOLS;
using namespace std;
using namespace EXTRAXS;

W_Decay_Virtual::W_Decay_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
  ME2_Base(pi, flavs) {
  Data_Reader reader(" ",";","#","=");
  // 0: Only QED plus soft real, 1: weak + QED + soft real, 2: weak (no soft real)
  double ew_prelim = reader.GetValue<int>("YFS_EW_CORRECTIONS",0);
  m_ew_scheme = reader.GetValue<int>("YFS_EW_SCHEME",2); 
  m_qed = 1;
  // Check born if necessary
  m_lo = reader.GetValue<int>("YFS_EXTRAXS_LO",0); 
  if (ew_prelim == 2) m_qed = 0;
  if (ew_prelim > 0) {
    m_ew = 1;
  }
  else m_ew = 0;
  if (m_ew == 0 && m_qed == 0) msg_Out() << "Turned every correction off. This should not happen!\n";
  msg_Debugging() << m_ew << " , " << m_qed << " , " << m_lo << " , " << ew_prelim;
  widthscheme=reader.GetValue<string>("WIDTH_SCHEME","CMS");
  // Infrared cutoff as used in WZGRAD
  // Conversion: omega_YFS = 0.5*deltas*shat
  m_deltas = reader.GetValue<double>("YFS_DELTAS",0.001);
  Complex I = Complex(0.,1.);
  // Save quark-flavours for determination of CKM element
  m_quark_flavs[0] = flavs[0];
  m_quark_flavs[1] = flavs[1];

  m_quark_masses[0] = (Flavour(kf_d).Mass()==0.?0.06984:Flavour(kf_d).Mass());
  m_quark_masses[1] = (Flavour(kf_u).Mass()==0.?0.06983:Flavour(kf_u).Mass());
  m_quark_masses[2] = (Flavour(kf_s).Mass()==0.?0.15:Flavour(kf_s).Mass());
  m_quark_masses[3] = (Flavour(kf_c).Mass()==0.?1.2:Flavour(kf_c).Mass());
  m_quark_masses[4] = (Flavour(kf_b).Mass()==0.?4.6:Flavour(kf_b).Mass());
  

  // Electroweak parameters, gauge boson masses
  double  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  double  MH  = Flavour(kf_h0).Mass();
  double  mt  = Flavour(kf_t).Mass();
  MW2 = pow(MW,2.);
  MZ2 = pow(MZ,2.);
  MH2 = pow(MH,2.);
  mt2 = pow(mt,2.);
  double  GH  = Flavour(kf_h0).Width();
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  GZ2 = GZ*GZ;
  GW2 = GW*GW;
  // if (MODEL::s_model->ScalarNumber("WidthScheme")) {
  muW2 = MW*(MW-I*GW);
  muZ2 = MZ*(MZ-I*GZ);
  muH2 = MH*(MH-I*GH);
  mut2 = mt*(mt-I*Flavour(kf_t).Width());
  if (widthscheme == "CMS") {
    m_cW2=muW2/muZ2;
    m_sW2=1.-m_cW2;
  }
  else if (widthscheme == "Fixed") {
    //muW2 = MW*MW*Complex(1.,0.);
    //muZ2 = MZ*MZ*Complex(1.,0.);
    //muH2 = MH*MH*Complex(1.,0.);
    //mut2 = mt*mt*Complex(1.,0.);
    m_cW2=MW2/MZ2;
    m_sW2=1.-m_cW2;
  }
  m_sW = sqrt(m_sW2);
  m_cW = sqrt(m_cW2);
  m_cW4 = m_cW2*m_cW2;
  m_sW4 = m_sW2*m_sW2;
  m_e = sqrt(4.*M_PI*AlphaQED()*CouplingFactor(0,1));
  m_alpha = AlphaQED()*CouplingFactor(0,1);

  // flavs is always ordered such that the neutrino sits at flavs[2],
  // the lepton at flavs[3]. Lepton charge determines the W charge
  // most easily. In decay setup (m_flavs), need to order such that
  // lepton comes before the neutrino though.
  m_flavs[0]  = (flavs[3].IsAnti()?Flavour(kf_Wplus,0):Flavour(kf_Wplus,1));
  m_flavs[1] = flavs[3]; 
  m_flavs[2] = flavs[2]; 
  m_cL = Complex(1.,0.);
  m_cR = Complex(0.,0.);


  Flavour IsoPartner;
  if (m_flavs[1].IsUptype()) {
    IsoPartner = ATOOLS::Flavour(m_flavs[1].Kfcode()-1,m_flavs[1].IsAnti());
  }
  else if (m_flavs[1].IsDowntype()) {
    IsoPartner = ATOOLS::Flavour(m_flavs[1].Kfcode()+1,m_flavs[1].IsAnti());
  }
  else {
    IsoPartner = ATOOLS::Flavour(0);
  }
  m_m2p = pow(IsoPartner.Mass(),2.);
  m_m2 = pow(m_flavs[1].Mass(),2.);
  m_Qf = m_flavs[1].Charge();
  m_If = m_flavs[1].IsoWeak();
  m_Qfp = IsoPartner.Charge();
  m_Ifp = IsoPartner.IsoWeak();

  msg->SetPrecision(16); 
  One = DivArrC(0.,0.,0.,1.,0.,0.);
  Zero = DivArrC(0.,0.,0.,0.,0.,0.);
}

W_Decay_Virtual::~W_Decay_Virtual() {
}

double W_Decay_Virtual::operator()(const Vec4D_Vector& momenta) {
  // Boost into CMS
  Poincare CMS = Poincare(momenta[0]+momenta[1]);
  Vec4D temp;
  temp = momenta[0]+momenta[1];
  CMS.Boost(temp);
  m_moms[0] = temp;
  temp = momenta[2];
  CMS.Boost(temp);
  m_moms[2] = temp;
  temp = momenta[3];
  CMS.Boost(temp);
  m_moms[1] = temp;
  m_masses[0] = m_moms[0].Abs();
  m_masses[1] = m_flavs[1].Mass();
  m_masses[2] = m_flavs[2].Mass();
  double sum = 0.;
  double s(0.),t(0.);
  s=(momenta[0]+momenta[1]).Abs2();
  // t is the invariant mass between incoming u-type quark and outgoing lepton
  // ordering of particles is odd in this respect hence this distinction
  if (m_quark_flavs[0].IsUptype() && m_flavs[2].Charge() != 0.) {
    t=(momenta[0]-momenta[2]).Abs2();
  }
  else if (m_quark_flavs[0].IsUptype() && m_flavs[2].Charge() == 0.) {
    t=(momenta[0]-momenta[3]).Abs2();
  }
  else if (m_quark_flavs[0].IsDowntype() && m_flavs[2].Charge() != 0.) {
    t=(momenta[1]-momenta[2]).Abs2();
  }
  else if (m_quark_flavs[0].IsDowntype() && m_flavs[2].Charge() == 0.) {
    t=(momenta[1]-momenta[3]).Abs2();
  }
  double CKM;
  // Determine relevant CKM element
  if ((m_quark_flavs[0].Kfcode() == kf_d && m_quark_flavs[1].Kfcode() == kf_u) ||
      (m_quark_flavs[0].Kfcode() == kf_u && m_quark_flavs[1].Kfcode() == kf_d)) CKM = 0.975;
  else if ((m_quark_flavs[0].Kfcode() == kf_s && m_quark_flavs[1].Kfcode() == kf_u) || 
	   (m_quark_flavs[0].Kfcode() == kf_u && m_quark_flavs[1].Kfcode() == kf_s)) CKM = 0.222;
  else if ((m_quark_flavs[0].Kfcode() == kf_d && m_quark_flavs[1].Kfcode() == kf_c) || 
	   (m_quark_flavs[0].Kfcode() == kf_c && m_quark_flavs[1].Kfcode() == kf_d)) CKM = 0.222;
  else if ((m_quark_flavs[0].Kfcode() == kf_s && m_quark_flavs[1].Kfcode() == kf_c) || 
	   (m_quark_flavs[0].Kfcode() == kf_c && m_quark_flavs[1].Kfcode() == kf_s)) CKM = 0.975;
  chi = sqr(t)/(sqr(s-MW2) + GW2*MW2);
  born = sqr(CKM)/(4.*abs(m_sW2)*abs(m_sW2))*chi*(1.-sqr(m_masses[1])/(2.*MW2)-pow(m_masses[1],4.)/(2.*sqr(MW2)));
  // Check born if necessary
  if (m_lo == 1) {
    m_res.Finite()=born;
    return sqr(4.*M_PI*m_alpha*CouplingFactor(0,1))*1./3.*born;
  }
  // Assemble virtual contribution
  DivArrC real = (m_qed?m_alpha*CouplingFactor(0,1)*YFS_Form_Factor()*One:Zero);
  DivArrC virt = (GetBeta_0_1()/GetBeta_0_0()+real)*born;
  // 1/epsIR
  m_res.IR()=virt.IR().real();
  // 1/epsIR2
  m_res.IR2()=virt.IR2().real();
  // finite
  m_res.Finite()=virt.Finite().real();
  return sqr(4.*M_PI*m_alpha*CouplingFactor(0,1))*1./3.*m_res.Finite().real();
}


Complex W_Decay_Virtual::InfraredSubtractedME_0_0() {
  // born ME
  Vec4C epsW = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_spins[1],epsW,2,m_spins[2],m_cR,m_cL);
}



DivArrC W_Decay_Virtual::InfraredSubtractedME_0_1() {
  // Setting up proper values for EW parameters:
  // Qf = 1; corresponding to If = 1/2 and m2 = ml2. Ifp = -1/2 and Qfp = 0, mfp = 0.
  // Need to offset the extra minus signs in .Charge() and .IsoWeak() for anti-particles
  if (!m_flavs[1].IsAnti()) {
    m_p1 = m_moms[1];
    m_p2 = m_moms[2];
    m_m2 = pow(m_masses[1],2.);
    m_m2p = pow(m_masses[2],2.);
    m_Qf = -m_flavs[1].Charge();
    m_Qfp = m_flavs[2].Charge();
    m_If = -m_flavs[1].IsoWeak();
    m_Ifp = m_flavs[2].IsoWeak();
  }
  else if (m_flavs[1].IsAnti()) {
    m_p1 = m_moms[1];
    m_p2 = m_moms[2];
    m_m2 = pow(m_masses[1],2.);
    m_m2p = pow(m_masses[2],2.);
    m_Qf = m_flavs[1].Charge();
    m_Qfp = -m_flavs[2].Charge();
    m_If = m_flavs[1].IsoWeak();
    m_Ifp = -m_flavs[2].IsoWeak();
  }
  Vec4C epsV = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  double s((m_moms[1]+m_moms[2]).Abs2());
  double mu2(s);
  m_mu2 = s;
  // Set up EW corrections as corrections to cR,cL
  // If only QED corrections, all weak form factors and QED contributions in CTs are neglected
  DivArrC term(0.,0.,0.,0.,0.,0.);

  Complex gfm((m_If-m_sW2*m_Qf)/(m_sW*m_cW)),gfpm((m_Ifp-m_sW2*m_Qfp)/(m_sW*m_cW));
  Complex pref(Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW));

  DivArrC cRDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(pref*(CT_R()));
  DivArrC cLDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(pref*(
							     CT_L()
							     +(m_Qf==-1?-1.:1.)*1./2.*FAn(s)
							     +((m_ew>=1)?(1./(4.*m_sW2*m_cW2)*(m_Ifp-m_sW2*m_Qfp)*(m_If-m_sW2*m_Qf)*FZa(s)
									  +(m_Qf==-1?1.:-1.)*1./m_sW2*(m_Ifp-m_sW2*m_Qfp-m_If+m_sW2*m_Qf)*FZn(s)):DivArrC(0.,0.,0.,0.,0.,0.))
							     ));
  DivArrC cRDivArrY = Zero;
  DivArrC cLDivArrY = Zero;
  // calculate infrared factor B
  DivArrC B(0.,0.,0.,0.,0.,0.);
  // m2*C_0(m2,m2,0,0,m2,m2) contains an IR divergence necessary to cancel all divergences even 
  // when m2 = 0 (C_0 ~ 1/m2)
  if (m_qed) {
    if (m_m2 != 0) {
      B += m_alpha*CouplingFactor(0,1)/M_PI*pref*(-0.5*(s-m_m2p+m_m2)*C_0(m_m2,m_m2p,s,0.,m_m2,muW2,mu2)
			      +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
			      +0.5*m_m2*C_0(m_m2,m_m2,0.,0.,m_m2,m_m2,mu2)
			      +0.25*B_0(m_m2p,s,m_m2,mu2)
			      -0.125*B_0(0.,m_m2,m_m2,mu2)
			      -0.125*B_0(0.,muW2,muW2,mu2))*
	XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
    }
    else if (m_m2 == 0) {
      B += m_alpha*CouplingFactor(0,1)/M_PI*pref*(-0.5*(s-m_m2p+m_m2)*C_0(m_m2,m_m2p,s,0.,m_m2,s,mu2)
			      +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
			      +0.25*DivArrC(0.,1.,0.,0.,0.,0.)
			      +0.25*B_0(m_m2p,s,m_m2,mu2)
			      -0.125*B_0(0.,m_m2,m_m2,mu2)
			      -0.125*B_0(0.,s,s,mu2))*
	XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
    }      
  }
  term = 
    DivArrC(XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.UV(),cLDivArr.UV()),
	    XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.IR(),cLDivArr.IR()),
	    XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.IR2(),cLDivArr.IR2()),
	    XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.Finite(),cLDivArr.Finite()),
	    XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.Epsilon(),cLDivArr.Epsilon()),
	    XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArr.Epsilon2(),cLDivArr.Epsilon2()))
    +DivArrC(XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.UV(),cLDivArrY.UV()),
	     XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.IR(),cLDivArrY.IR()),
	     XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.IR2(),cLDivArrY.IR2()),
	     XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Finite(),cLDivArrY.Finite()),
	     XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Epsilon(),cLDivArrY.Epsilon()),
	     XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Epsilon2(),cLDivArrY.Epsilon2()))
    +B
    ; 
  return term;
}

Complex W_Decay_Virtual::GetBeta_0_0() {
  // Assemble born level matrix element squared
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
        m_spins[0] = k;
        m_spins[1] = j;
        m_spins[2] = i;
        Complex M_0_0 = InfraredSubtractedME_0_0();
        sum = sum + (M_0_0*conj(M_0_0))// .real()
	  ;
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  return sum;
}

DivArrC W_Decay_Virtual::GetBeta_0_1() {
  // Assemble virtual matrix element squared
  DivArrC sum(0.,0.,0.,0.,0.,0.);
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
	m_spins[0] = k;
	m_spins[1] = j;
	m_spins[2] = i;
	Complex M_0_0 = InfraredSubtractedME_0_0();
	DivArrC M_0_1 = InfraredSubtractedME_0_1();
	sum = sum + (M_0_0*conj(M_0_1)+conj(M_0_0)*M_0_1);
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  return sum;
}


// EW form factors (number relate to equations in Bardin:1999ak)
// WWy vertex (5.606)
DivArrC W_Decay_Virtual::FAn(const double& s) {
  Complex w(-s/muW2);
  return m_Qf*(s*C_0(m_m2,m_m2p,s,0.,m_m2,muW2,m_mu2)+B_0(m_m2,m_m2,0.,m_mu2))
    -m_Qfp*(s*C_0(m_m2,m_m2p,s,muW2,m_m2p,0.,m_mu2)+B_0(m_m2p,m_m2p,0.,m_mu2))
    +(m_Qf-m_Qfp)/2.*(-2.*B_0(s,muW2,0.,m_mu2)+3.*B_0(0.,0.,muW2,m_mu2));
}

// abelian Z-exchange (5.612 and 5.570)
DivArrC W_Decay_Virtual::FZa(const double& s) { 
  Complex z(-s/muZ2);
  return 2.*muZ2/z*pow((1.-z),2.)*C_0(0.,0.,s,0.,muZ2,0.,m_mu2)
    + B_0(s,0.,0.,m_mu2)
    + (2./z-4.)*(B_0(s,0.,0.,m_mu2)-B_0(0.,0.,muZ2,m_mu2))
    -2.*DivArrC(0.,0.,0.,1.,0.,0.);
}

// WWZ vertex (5.615)
DivArrC W_Decay_Virtual::FZn(const double& s) { 
  Complex w(-s/muW2), z(-s/muZ2);
  return 0.5*(-((1./w-1.)/m_cW2-1.)*muW2*C_0(0.,0.,s,muW2,0.,muZ2,m_mu2)
	      +0.5*(1./z+1./w-1.)*B_0(s,muW2,muZ2,m_mu2)
	      -(0.5/z-1.)/muZ2*A_0(muZ2,m_mu2)
	      -(0.5/w-1.)/muW2*A_0(muW2,m_mu2)); /// check signs
}

// Left handed counterterm
DivArrC W_Decay_Virtual::CT_L() { 
  if (m_ew == 0) return Complex(1.,0.)*(dZe() + 0.5*dZW() + 0.5*(dZfermL(m_m2p,m_m2,m_Qfp,m_Ifp) + dZfermL(m_m2,m_m2p,m_Qf,m_If)));
  return dZe() + m_cW2/m_sW2*dcw() + 0.5*dZW() + 0.5*(dZfermL(m_m2p,m_m2,m_Qfp,m_Ifp) + dZfermL(m_m2,m_m2p,m_Qf,m_If));
}

DivArrC W_Decay_Virtual::CT_R() { 
  // No right handed coupling at LO, hence no counterterm contribution
  return Zero;
}


// YFS form factor calculation
double W_Decay_Virtual::YFS_Form_Factor() {
  double s = (m_moms[1]+m_moms[2]).Abs2();
  double mu2 = s;
  m_omega = m_deltas*sqrt(s)/2.;
  m_p1 = m_moms[0];
  m_p2 = m_moms[1];
  m_m_1 = m_moms[0].Abs();//m_moms[1].Abs();
  m_m_2 = m_flavs[1].Mass();//m_moms[2].Abs();
  // PRINT_VAR("\n" << m_p1 << "\n" << m_p2 << "\n" << m_m_1 << "\n" << m_m_2);
  // if (m_moms[0][0] >= m_moms[1][0]) {
  //   std::swap(m_p1,m_p2);
  //   std::swap(m_m_1,m_m_2);
  // }
  // m_x1  = 0.;
  // m_x2  = 0.;
  // return 1./M_PI*(log(m_p1[0]*m_p2[0]/sqr(m_omega)) + 0.5*(m_p1*m_p2)*IntP1()
  // 		  - 0.5*(m_p1*m_p2)*IntE() + 0.25*IntP2() + G(1.) + G(-1.)
  // 		  - (m_p1*m_p2)*IntG());



  /// Energy cutoff is applied in the partonic centre of mass frame, i. e. the rest frame of the 
  /// W boson. This is not the same as in PHOTONS, hence need different expression for YFS form factor
  /// Use expression at the end of hep-ph/0302065
  return 1./M_PI*(2.*(log(m_m_1/m_m_2)-1.)*log(2.*m_omega/m_m_1) + 0.5*log(m_m_1/m_m_2) - 0.5 - sqr(M_PI)/6.);
}

double W_Decay_Virtual::IntP1() {
  return 0.;
}

double W_Decay_Virtual::IntE() {
  // E1 = E2
  // PRINT_VAR("\n" << 1./(m_x1-m_x2)*log((m_p1[0]+m_p2[0])/(2*m_omega))
  // 	    *log(abs(((1.-m_x1)*(1.+m_x2))/((1.+m_x1)*(1.-m_x2)))) << "\n" << -0.5*log(2.*(m_p1*m_p2)/sqr(m_m_1))*log(m_p1[0]*m_p2[0]/sqr(m_omega)));
  if (abs(m_p1[0]-m_p2[0]) > 1E-6) {
    double xE = - (m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
    double xp = - (m_m_1*m_m_1+m_m_2*m_m_2)/(m_m_1*m_m_1-m_m_2*m_m_2);
    // xE = xp
    if (abs(xE-xp) < 1E-6)
      return 4./(m_m_2*m_m_2-m_m_1*m_m_1)
	*(log((m_p2[0]-m_p1[0])/(2.*m_omega))*log(abs((xp+1.)/(xp-1.)))
	  - (1./2.)*(sqr(log(xp-1.))-sqr(log(xp+1.))));
    // xE > xp
    else if (xE > xp)
      return 4./(m_m_2*m_m_2-m_m_1*m_m_1)
	*(log((m_p2[0]-m_p1[0])/(2.*m_omega))*log(abs((xp+1.)/(xp-1.)))
	  + log(xE-xp)*log(abs((xp+1.)/(xp-1.)))
	  + DiLog((xp-1.)/(xp-xE)) - DiLog((xp+1.)/(xp-xE)));
    // xE < xp
    else if (xp > xE) {
      double xi = (m_p1[0]-m_p2[0])/(2*m_p2[0]);
      double yp = 1.+xi*(1.+xp);
      return 4./(m_m_2*m_m_2-m_m_1*m_m_1)
	*((log(m_p2[0]/m_omega)+log(abs(yp)))*log(abs((xp+1.)/(xp-1.)))
	  - (1./2.)*sqr(log(abs(yp/(yp-1.-2.*xi))))
	  + (1./2.)*sqr(log(abs(yp/(yp-1.))))
	  + DiLog(yp/(yp-1.)) - DiLog(yp/(yp-1.-2.*xi)));
    }
    else {
      msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
      return 0.;
    }
  }
  else {
    msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
    return 0.;
  }
}

// int dx ln (px'²/m1m2)                                                  (C.55)
double W_Decay_Virtual::IntP2() {
  if (abs(m_m_1*m_m_1 - m_m_2*m_m_2) > 1E-6) {
    double xp = -(m_m_1*m_m_1+m_m_2*m_m_2)/(m_m_1*m_m_1 - m_m_2*m_m_2);
    return 2.*log(abs(m_m_1*m_m_1-m_m_2*m_m_2)/(2.*m_m_1*m_m_2))
      + log(abs(1.-xp*xp)) + xp*log(abs((1.+xp)/(1.-xp))) - 2.;
  }
}

// G(x)
double W_Decay_Virtual::G(double x) {
  Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
  double b  = CalculateBeta(px);
  double r(0.);
  if (b == 0.)      r = 1.-LOG_2;
  else if (b == 1.) r = 0.;
  else              r = ((1.-b)/(2.*b)*log((1.+b)/(1.-b))+log((1.+b)/2.));
  return r;
}

// int dx G(x)/px²                                                        (C.88)
double W_Decay_Virtual::IntG() {
  // if dipole in its CMS
  if ((abs((m_p1-m_p2).Abs2()) < 1E-6) &&
      (m_p1.Abs2()/m_p2.Abs2() < 1E-3)) {
    return 2./m_p2.Abs2()*(3./12.*M_PI*M_PI+DiLog(-2.));
  }  
  else {
    msg_Out() << METHOD << "whut?\n";
  }
}

double W_Decay_Virtual::CalculateBeta(const Vec4D& p) {
  return (Vec3D(p).Abs()/p[0]);
}





// Ingredients to set up counterterms. Taken from Denner:1991kt
DivArrC W_Decay_Virtual::dZfermL(const double& m2, const double& m2p,
				 const double& Qf, const double& If) {
  // Fermion wavefunction counterterm, contribution to PL = (1 - gamma5)/2
  if (m2 == 0.) {
    return -(Sigma_ferm_L(m2,m2,m2p,Qf,If));
  }
  else 
    return -(Sigma_ferm_L(m2,m2,m2p,Qf,If))
      -(m2*(dSigma_ferm_R(m2,m2,m2p,Qf,If)
	    +dSigma_ferm_L(m2,m2,m2p,Qf,If)
	    +2.*dSigma_ferm_S(m2,m2,m2p,Qf,If)));
}

DivArrC W_Decay_Virtual::dZfermR(const double& m2, const double& m2p,
				 const double& Qf, const double& If) {
  // Fermion wavefunction counterterm, contribution to PR = (1 + gamma5)/2
  if (m2 == 0.) {
    return -(Sigma_ferm_R(m2,m2,m2p,Qf,If));
  }
  else 
    return -(Sigma_ferm_R(m2,m2,m2p,Qf,If))
      -m2*(dSigma_ferm_R(m2,m2,m2p,Qf,If)
	   +dSigma_ferm_L(m2,m2,m2p,Qf,If)
	   +2.*dSigma_ferm_S(m2,m2,m2p,Qf,If));
}

DivArrC W_Decay_Virtual::dm(const double& m2, const double& m2p,
			    const double& Qf, const double& If) {
  return sqrt(m2)/2.*(Sigma_ferm_L(m2,m2,m2p,Qf,If) 
		      + Sigma_ferm_R(m2,m2,m2p,Qf,If) 
		      + 2.*Sigma_ferm_S(m2,m2,m2p,Qf,If));
}

// Fermion self energies
DivArrC W_Decay_Virtual::Sigma_ferm_R(const double& p2,
				      const double& m2, const double& m2p,
				      const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf);
  return (m_qed?-1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One):Zero)
    +(m_ew?(-1./4.*pow(gplus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	    -1./(16.*m_sW2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2/muW2*B_1(p2,m2p,muW2,m_mu2)):Zero);
}

DivArrC W_Decay_Virtual::Sigma_ferm_L(const double& p2,
				      const double& m2, const double& m2p,
				      const double& Qf, const double& If)
{
  Complex gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One):Zero)
    +(m_ew?(-1./4.*pow(gminus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	    -1./(16.*m_sW2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*((2.+m2p/muW2)*B_1(p2,m2p,muW2,m_mu2)+1.*One)):Zero);
}


DivArrC W_Decay_Virtual::Sigma_ferm_S(const double& p2,
				      const double& m2, const double& m2p,
				      const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf),gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*(4.*B_0(p2,m2,0.,m_mu2)-2.*One):Zero)
    +(m_ew?(-1./4.*gplus*gminus*(4.*B_0(p2,m2,muZ2,m_mu2)-2.*One)
	    -1./(16.*m_sW2)*m2/muW2*(B_0(p2,m2,muZ2,m_mu2)-B_0(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2p/muW2*B_0(p2,m2p,muW2,m_mu2)):Zero);
}

// Derivatives of fermion self energies
DivArrC W_Decay_Virtual::dSigma_ferm_R(const double& p2,
				       const double& m2, const double& m2p,
				       const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf);
  return (m_qed?-1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2):Zero)
    +(m_ew?(-1./4.*pow(gplus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sW2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2/muW2*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}

DivArrC W_Decay_Virtual::dSigma_ferm_L(const double& p2,
				       const double& m2, const double& m2p,
				       const double& Qf, const double& If)
{
  Complex gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2):Zero)
    +(m_ew?(-1./4.*pow(gminus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sW2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*(2.+m2p/muW2)*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}


DivArrC W_Decay_Virtual::dSigma_ferm_S(const double& p2,
				       const double& m2, const double& m2p,
				       const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf),gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*4.*B_0p(p2,m2,0.,m_mu2):Zero)
    +(m_ew?(-1./4.*gplus*gminus*4.*B_0p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sW2)*m2/muW2*(B_0p(p2,m2,muZ2,m_mu2)-B_0p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2p/muW2*B_0p(p2,m2p,muW2,m_mu2)):Zero);
}


// W Boson self energy
DivArrC W_Decay_Virtual::Sigma_W(const double& p2)
{
  // if p2 == 0 pieces 1/p2*(B0(p2,...)-B0(0,...)) become dB0(0,...)
  
  int kfl[] = {11,12,13,14,15,16}; // lepton IDs
  int kfq[] = {1,2,3,4,5,6}; // quark IDs
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  for (int i = 0; i < 3; ++i) {
    Flavour flav(Flavour(kfl[2*i]));
    double m2(pow(flav.Mass(),2.));
    fermDenner += 1./(12.*m_sW2)*(-(p2-m2/2.)*B_0(p2,0.,m2,m_mu2) + 1./3.*p2*One
				  +m2*B_0(0.,m2,m2,m_mu2) 
				  + (p2!=0.?pow(m2,2.)/(2.*p2)*(B_0(p2,0.,m2,m_mu2) - B_0(0.,0.,m2,m_mu2)):pow(m2,2.)/2.*B_0p(p2,0.,m2,m_mu2)));
  }
  for (int i = 0; i < 3; ++i) {
    Flavour flavd(Flavour(kfq[2*i])),flavu(Flavour(kfq[2*i+1]));
    // make sure to use correct value of mt2 for comparison
    double m2d(sqr(m_quark_masses[2*i])), m2u(((i==2)?mt2:sqr(m_quark_masses[2*i+1])));//double m2d(pow(flavd.Mass(),2.)),m2u(((i==2)?mt2:pow(flavu.Mass(),2.)));
    fermDenner += 1./(4.*m_sW2)*(-(p2-(m2d+m2u)/2.)*B_0(p2,m2u,m2d,m_mu2)+1./3.*p2*One
				 +m2u*B_0(0.,m2u,m2u,m_mu2)+m2d*B_0(0.,m2d,m2d,m_mu2)
				 +(p2!=0.?pow(m2u-m2d,2.)/(2.*p2)*(B_0(p2,m2u,m2d,m_mu2)-B_0(0.,m2u,m2d,m_mu2)):pow(m2u-m2d,2.)/2.*B_0p(p2,m2u,m2d,m_mu2)));
  }
  return -(m_ew?fermDenner:Zero)
    - (1./6.*(((2.*muW2+5.*p2)*B_0(p2,muW2,0.,m_mu2) - 2.*muW2*B_0(0.,muW2,muW2,m_mu2)
	      -(p2!=0.?
		(pow(muW2,2.)/p2*(B_0(p2,muW2,0.,m_mu2)
				  -B_0(0.,muW2,0.,m_mu2)))
		:(pow(muW2,2.)*B_0p(p2,muW2,0.,m_mu2)))
		     +1./3.*p2*One))
       +(m_ew?(1./(48.*m_sW2)*(((40.*m_cW2-1.)*p2
				+(16.*m_cW2+54.-10./m_cW2)*muW2)*B_0(p2,muW2,muZ2,m_mu2)
			       -(16.*m_cW2+2.)*(muW2*B_0(0.,muW2,muW2,m_mu2)
						+muZ2*B_0(0.,muZ2,muZ2,m_mu2))
			       +(4.*m_cW2-1.)*2./3.*p2*One
			       -(8.*m_cW2+1.)*(p2!=0.?
					       (pow(muW2-muZ2,2.)/p2*(B_0(p2,muW2,muZ2,m_mu2)
								      -B_0(0.,muW2,muZ2,m_mu2))):
					       (pow(muW2-muZ2,2.)*B_0p(p2,muW2,muZ2,m_mu2)))
			       )
	       +1./(48.*m_sW2)*((2.*muH2-10.*muW2-p2)*B_0(p2,muW2,muH2,m_mu2)
				-2.*muW2*B_0(0.,muW2,muW2,m_mu2)
				-2.*muH2*B_0(0.,muH2,muH2,m_mu2)
				-(p2!=0.?
				  (pow(muW2-muH2,2.)/p2*(B_0(p2,muW2,muH2,m_mu2)
							 -B_0(0.,muW2,muH2,m_mu2))):
				  pow(muW2-muH2,2.)*B_0p(p2,muW2,muH2,m_mu2))
				-2./3.*p2*One)):Zero));
}


DivArrC W_Decay_Virtual::dSigma_W(const double& p2)
{
  //Denner
  int kfl[] = {11,12,13,14,15,16}; // lepton IDs
  int kfq[] = {1,2,3,4,5,6}; // quark IDs
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  for (int i = 0; i < 3; ++i) {
    Flavour flav(Flavour(kfl[2*i]));
    double m2(pow(flav.Mass(),2.));
    fermDenner += 1./(12.*m_sW2)*(-(p2-m2/2.)*B_0p(p2,0.,m2,m_mu2)
				  -B_0(p2,0.,m2,m_mu2)
				  + 1./3.*One
				  + pow(m2,2.)/(2.*p2)*B_0p(p2,0.,m2,m_mu2)
				  - pow(m2,2.)/(2.*pow(p2,2.))*(B_0(p2,0.,m2,m_mu2) - B_0(0.,0.,m2,m_mu2)));
  }
  for (int i = 0; i < 3; ++i) {
    Flavour flavd(Flavour(kfq[2*i])),flavu(Flavour(kfq[2*i+1]));
    // make sure to use correct value of mt2 for comparison
    double m2d(sqr(m_quark_masses[2*i])), m2u(((i==2)?mt2:sqr(m_quark_masses[2*i+1])));//double m2d(pow(flavd.Mass(),2.)),m2u(((i==2)?mt2:pow(flavu.Mass(),2.)));
    fermDenner += 1./(4.*m_sW2)*(-(p2-(m2d+m2u)/2.)*B_0p(p2,m2u,m2d,m_mu2)
				 -B_0(p2,m2u,m2d,m_mu2)
				 +1./3.*One
				 +pow(m2u-m2d,2.)/(2.*p2)*B_0p(p2,m2u,m2d,m_mu2)
				 -pow(m2u-m2d,2.)/(2.*pow(p2,2.))*(B_0(p2,m2u,m2d,m_mu2)-B_0(0.,m2u,m2d,m_mu2)));
  }
  return -(m_ew?fermDenner:Zero)
    - (1./6.*(((2.*muW2+5.*p2)*B_0p(p2,muW2,0.,m_mu2)
	      +5.*B_0(p2,muW2,0.,m_mu2)
	      -pow(muW2,2.)/p2*B_0p(p2,muW2,0.,m_mu2)
	      +pow(muW2/p2,2.)*(B_0(p2,muW2,0.,m_mu2)
				-B_0(0.,muW2,0.,m_mu2))
		     +1./3.*One))
       +(m_ew?(1./(48.*m_sW2)*(((40.*m_cW2-1.)*p2
				+(16.*m_cW2+54.-10./m_cW2)*muW2)*B_0p(p2,muW2,muZ2,m_mu2)
			       +(40.*m_cW2-1.)*B_0(p2,muW2,muZ2,m_mu2)
			       +(4.*m_cW2-1.)*2./3.*One
			       -(8.*m_cW2+1.)*pow(muW2-muZ2,2.)/p2*B_0p(p2,muW2,muZ2,m_mu2)
			       +(8.*m_cW2+1.)*pow((muW2-muZ2)/p2,2.)*(B_0(p2,muW2,muZ2,m_mu2)
								      -B_0(0.,muW2,muZ2,m_mu2)))
	       +1./(48.*m_sW2)*((2.*muH2-10.*muW2-p2)*B_0p(p2,muW2,muH2,m_mu2)
				-B_0(p2,muW2,muH2,m_mu2)
				-pow(muW2-muH2,2.)/p2*B_0p(p2,muW2,muH2,m_mu2)
				+pow((muW2-muH2)/p2,2.)*(B_0(p2,muW2,muH2,m_mu2)
							 -B_0(0.,muW2,muH2,m_mu2))
				-2./3.*One)):Zero));
}



// Z Boson self energy
DivArrC W_Decay_Virtual::Sigma_ZZ(const double& p2)
{
  if (m_ew == 0) return Zero; 
  int N_f(12);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  // Denner
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2;
    if (i <= 5) m2 = (i==5?mt2:sqr(m_quark_masses[i]));
    else m2 = pow(flav.Mass(),2.);
    Complex gplus(-m_sW/m_cW*flav.Charge()),gminus((flav.IsoWeak()-m_sW2*flav.Charge())/(m_sW*m_cW));
    fermDenner += 1./6.*((flav.IsQuark())?3.:1.)
      *((pow(gplus,2.)+pow(gminus,2.))*(-(p2+2.*m2)*B_0(p2,m2,m2,m_mu2) 
					+ 2.*m2*B_0(0.,m2,m2,m_mu2) 
					+ 1./3.*p2*One)
	+3./(4.*m_sW2*m_cW2)*m2*B_0(p2,m2,m2,m_mu2));
  }
  return -fermDenner 
    -(1./(24.*m_sW2*m_cW2)*(((18.*m_cW4+2.*m_cW2-1./2.)*p2
			     +(24.*m_cW4+16.*m_cW2-10.)*muW2)*B_0(p2,muW2,muW2,m_mu2)
			    -(24.*m_cW4-8.*m_cW2+2.)*muW2*B_0(0.,muW2,muW2,m_mu2)
			    +(4.*m_cW2-1.)*1./3.*p2*One)
      +1./(48.*m_sW2*m_cW2)*((2.*muH2-10.*muZ2-p2)*B_0(p2,muZ2,muH2,m_mu2)
			     -2.*muZ2*B_0(0.,muZ2,muZ2,m_mu2) - 2.*muH2*B_0(0.,muH2,muH2,m_mu2)
			     -pow(muZ2-muH2,2.)/p2*(B_0(p2,muZ2,muH2,m_mu2)-B_0(0.,muZ2,muH2,m_mu2))
			     -2./3.*p2*One));
}

// Derivative of Z Boson self energy
DivArrC W_Decay_Virtual::dSigma_ZZ(const double& p2)
{
  if (m_ew == 0) return Zero; 
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; //fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2;
    if (i <= 5) m2 = (i==5?mt2:sqr(m_quark_masses[i]));
    else m2 = pow(flav.Mass(),2.);
    Complex gplus(-m_sW/m_cW*flav.Charge()),gminus((flav.IsoWeak()-m_sW2*flav.Charge())/(m_sW*m_cW));
    fermDenner += 1./6.*((flav.IsQuark())?3.:1.)
      *((pow(gplus,2.)+pow(gminus,2.))*(-(p2+2.*m2)*B_0p(p2,m2,m2,m_mu2) 
					-B_0(p2,m2,m2,m_mu2)
					+ 1./3.*One)
	+3./(4.*m_sW2*m_cW2)*m2*B_0p(p2,m2,m2,m_mu2));
  }
  return -fermDenner
    -(1./(24.*m_sW2*m_cW2)*(((18.*m_cW4+2.*m_cW2-1./2.)*p2
					+(24.*m_cW4+16.*m_cW2-10.)*muW2)*B_0p(p2,muW2,muW2,m_mu2)
				       +(18.*m_cW4+2.*m_cW2-1./2.)*B_0(p2,muW2,muW2,m_mu2)
				       +(4.*m_cW2-1.)*1./3.*One)
		 +1./(48.*m_sW2*m_cW2)*((2.*muH2-10.*muZ2-p2)*B_0p(p2,muZ2,muH2,m_mu2)
					-B_0(p2,muZ2,muH2,m_mu2)
					-pow(muZ2-muH2,2.)/p2*B_0p(p2,muZ2,muH2,m_mu2)
					+pow(muZ2-muH2,2.)/pow(p2,2.)*(B_0(p2,muZ2,muH2,m_mu2)-B_0(0.,muZ2,muH2,m_mu2))
					-2./3.*One));
}

// Higgs Boson self energy
DivArrC W_Decay_Virtual::Sigma_H(const double& p2)
{
  // Denner
  if (m_ew == 0) return Zero; 
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.); 
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2;
    if (i <= 5) m2 = (i==5?mt2:sqr(m_quark_masses[i]));
    else m2 = pow(flav.Mass(),2.);
    fermDenner += 1./(8.*m_sW2*muW2)*(flav.IsQuark()?3.:1.)*m2*(2.*A_0(m2,m_mu2) + (4.*m2-p2)*B_0(p2,m2,m2,m_mu2));
  }
  return -fermDenner 
    - (-1./(8.*m_sW2)*((6.*muW2-2.*p2+pow(muH2,2.)/(2.*muW2))*B_0(p2,muW2,muW2,m_mu2)
				  +(3.+muH2/(2.*muW2))*A_0(muW2,m_mu2) - 6.*muW2)
		  -1./(16.*m_sW2*m_cW2)*((6.*muZ2-2.*p2+pow(muH2,2.)/(2.*muZ2))*B_0(p2,muZ2,muZ2,m_mu2)
					 +(3.+muH2/(2.*muZ2))*A_0(muZ2,m_mu2) - 6.*muZ2)
		  -3./(32.*m_sW2)*(3.*pow(muH2,2.)/muW2*B_0(p2,muH2,muH2,m_mu2)
				   +muH2/muW2*A_0(muH2,m_mu2)));
}

// Derivative of Higgs Boson self energy
DivArrC W_Decay_Virtual::dSigma_H(const double& p2)
{
  // Denner
  if (m_ew == 0) return Zero; 
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2;
    if (i <= 5) m2 = (i==5?mt2:sqr(m_quark_masses[i]));
    else m2 = pow(flav.Mass(),2.);
    fermDenner += 1./(8.*m_sW2*muW2)*(flav.IsQuark()?3.:1.)*m2*((4.*m2-p2)*B_0p(p2,m2,m2,m_mu2)
							       -B_0(p2,m2,m2,m_mu2));
  }
  return -fermDenner 
    - (-1./(8.*m_sW2)*((6.*muW2-2.*p2+pow(muH2,2.)/(2.*muW2))*B_0p(p2,muW2,muW2,m_mu2)
				  -2.*B_0(p2,muW2,muW2,m_mu2))
		  -1./(16.*m_sW2*m_cW2)*((6.*muZ2-2.*p2+pow(muH2,2.)/(2.*muZ2))*B_0p(p2,muZ2,muZ2,m_mu2)
					 -2.*B_0(p2,muZ2,muZ2,m_mu2))
		  -3./(32.*m_sW2)*3.*pow(muH2,2.)/muW2*B_0p(p2,muH2,muH2,m_mu2));
}


// Derivative of Photon self energy
DivArrC W_Decay_Virtual::dSigma_GamGam(const double& p2)
{
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2;
    if (i <= 5) m2 = (i==5?mt2:sqr(m_quark_masses[i]));
    else m2 = pow(flav.Mass(),2.);
    if(m_ew) fermDenner += 1./3.*(flav.IsQuark()?3.:1.)*pow(flav.Charge(),2.)*(-B_0(p2,m2,m2,m_mu2)
									       -(p2+2.*m2)*B_0p(p2,m2,m2,m_mu2)
									       +1./3.*One);
    else fermDenner += 1./3.*(flav.IsQuark()?3.:1.)*pow(flav.Charge(),2.)*(-B_0(p2,m2,m2,m_mu2)
									   +1./3.*One);
  }
  return - fermDenner
    - (m_ew?(1./4.*((3.*p2+4.*muW2)*B_0p(p2,muW2,muW2,m_mu2) 
				   +3.*B_0(p2,muW2,muW2,m_mu2))):Zero);
}

// Derivative of photon self energy at p2 = 0., including some terms to resum leading logarithms
// from fermion contributions in alpha(0) scheme - see OpenLoops implementation
DivArrC W_Decay_Virtual::dSigma_GamGam0()
{
  DivArrC PiAALightZ = Zero;
  double m_u2 = pow(Flavour(1).Mass(),2.);
  double md2 = pow(Flavour(2).Mass(),2.);
  double ms2 = pow(Flavour(3).Mass(),2.);
  double mc2 = pow(Flavour(4).Mass(),2.);
  double mb2 = pow(Flavour(5).Mass(),2.);
  double me2 = pow(Flavour(11).Mass(),2.);
  double mm2 = pow(Flavour(13).Mass(),2.);
  double mtau2 = pow(Flavour(15).Mass(),2.);
  PiAALightZ += 1./3.*
    (2.*me2*B_0(0.,me2,me2,m_mu2) - (2.*me2+muZ2)*B_0(MZ2,me2,me2,m_mu2)
     +2.*mm2*B_0(0.,mm2,mm2,m_mu2) - (2.*mm2+muZ2)*B_0(MZ2,mm2,mm2,m_mu2)
     +2.*mtau2*B_0(0.,mtau2,mtau2,m_mu2) - (2.*mtau2+muZ2)*B_0(MZ2,mtau2,mtau2,m_mu2)
     +4./3.*(2.*m_u2*B_0(0.,m_u2,m_u2,m_mu2) - (2.*m_u2+muZ2)*B_0(MZ2,m_u2,m_u2,m_mu2))
     +1./3.*(2.*md2*B_0(0.,md2,md2,m_mu2) - (2.*md2+muZ2)*B_0(MZ2,md2,md2,m_mu2)) 
     +1./3.*(2.*ms2*B_0(0.,ms2,ms2,m_mu2) - (2.*ms2+muZ2)*B_0(MZ2,ms2,ms2,m_mu2)) 
     +4./3.*(2.*mc2*B_0(0.,mc2,mc2,m_mu2) - (2.*mc2+muZ2)*B_0(MZ2,mc2,mc2,m_mu2)) 
     +1./3.*(2.*mb2*B_0(0.,mb2,mb2,m_mu2) - (2.*mb2+muZ2)*B_0(MZ2,mb2,mb2,m_mu2)) 
     +muZ2*(1.*One+11./27.*One));
  DivArrC dSigmaTop = 4./9.*(1./3.*One-B_0(0.,mut2,mut2,m_mu2)+2.*mt2*B_0p(0.,mut2,mut2,m_mu2));
  double alphaQED_0 = 1./137.035999074;
  double alphaQED_MZ = m_alpha;
  double dAlphaQED_MZ = M_PI/(m_alpha*CouplingFactor(0,1))*(1.-alphaQED_0/alphaQED_MZ);
  return - (PiAALightZ/muZ2 + dSigmaTop)
    - (m_ew?(1./4.*(4.*muW2*B_0p(0.,muW2,muW2,m_mu2) 
  				   +3.*B_0(0.,muW2,muW2,m_u2))):Zero)
    +dAlphaQED_MZ*One;
}

// Z-Photon transition
DivArrC W_Decay_Virtual::Sigma_ZA(const double& p2)
{
  if (m_ew == 0) return Zero;
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:sqr(m_quark_masses[i]))),ch(flav.Charge());//double m2((i==5?mt2:pow(flav.Mass(),2.))),ch(flav.Charge());
    Complex gplus(-m_sW/m_cW*flav.Charge()),gminus((flav.IsoWeak()-m_sW2*flav.Charge())/(m_sW*m_cW));
    fermDenner += 1./6.*(flav.IsQuark()?3.:1.)*(-flav.Charge())*(gplus+gminus)*(-(p2+2.*m2)*B_0(p2,m2,m2,m_mu2)
								      +2.*m2*B_0(0.,m2,m2,m_mu2)
										+1./3.*p2*One);
  }
  return -fermDenner
    + (m_ew?(1./(12.*m_sW*m_cW)*(((9.*m_cW2+1./2.)*p2+(12.*m_cW2+4.)*muW2)*B_0(p2,muW2,muW2,m_mu2)
				      -(12.*m_cW2-2.)*muW2*B_0(0.,muW2,muW2,m_mu2)
				      +1./3.*p2*One)):Zero);
}

// Wavefunction Counterterms
DivArrC W_Decay_Virtual::dZW()
{
  return -dSigma_W(MW2);
}

DivArrC W_Decay_Virtual::dZZZ()
{
  return -dSigma_ZZ(MZ2);
}

DivArrC W_Decay_Virtual::dZAA()
{
  if (m_ew_scheme == 1 || m_ew_scheme == 3) return -dSigma_GamGam(MZ2);
  else if (m_ew_scheme == 2) return -dSigma_GamGam0();
}

DivArrC W_Decay_Virtual::dZAZ()
{
  if (widthscheme == "Fixed") return -Complex(1.,0.)*2.*real(Sigma_ZA(MZ2))/MZ2;
  else if (widthscheme == "CMS") return -2.*Sigma_ZA(MZ2)/MZ2 + (muZ2/MZ2 - 1.)*dZZA();
}

DivArrC W_Decay_Virtual::dZZA()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*2.*Sigma_ZA(0.)/muZ2;
  else if (widthscheme == "CMS") return 2.*Sigma_ZA(0.)/muZ2;
}

DivArrC W_Decay_Virtual::dZH()
{
  return -dSigma_H(MH2);
}

// Mass Counterterms
DivArrC W_Decay_Virtual::dMW2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_W(MW2));
  else if (widthscheme == "CMS") return Sigma_W(MW2) + (muW2 - MW2)*dSigma_W(MW2);
}

DivArrC W_Decay_Virtual::dMZ2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_ZZ(MZ2));
  else if (widthscheme == "CMS") return Sigma_ZZ(MZ2) + (muZ2-MZ2)*dSigma_ZZ(MZ2);
}

DivArrC W_Decay_Virtual::dMH2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_H(MH2));
  else if (widthscheme == "CMS") return Sigma_H(MH2) + (muH2-MH2)*dSigma_H(MH2);
}

// Counterterm for cos(thetaW)
DivArrC W_Decay_Virtual::dcw()
{
  if (widthscheme == "Fixed") return 1./2.*(dMW2()/MW2-dMZ2()/MZ2);
  else if (widthscheme == "CMS") return 1./2.*(dMW2()/muW2-dMZ2()/muZ2);
}

// Charge renormalization counterterm
DivArrC W_Decay_Virtual::dZe()
{
  switch (m_ew_scheme) {
  case 1: {
    if (widthscheme == "Fixed") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/MZ2;  /// This is correct for \alpha(MZ) scheme
    if (widthscheme == "CMS") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/muZ2;  /// This is correct for \alpha(MZ) scheme
    break;
  }
  case 2: {
    return 1./2.*dSigma_GamGam0() - m_sW/m_cW*Sigma_ZA(0.)/muZ2;
    break;
  }
  case 3: {
    return -m_cW2/m_sW2*dcw() + 1./2.*(Sigma_W(MW2)-Sigma_W(0.))/muW2 
      - 1./(m_cW*m_sW)*Sigma_ZA(0.)/muZ2 - One*1./8.*1./m_sW2*(6.+(7.-4.*m_sW2)/(2.*m_sW2)*log(m_cW2));
    break;
  }
  default: { 
    if (widthscheme == "Fixed") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/MZ2;  /// This is correct for \alpha(MZ) scheme
    if (widthscheme == "CMS") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/muZ2;  /// This is correct for \alpha(MZ) scheme
    break;
  }
  }
}


DECLARE_TREEME2_GETTER(W_Decay_Virtual,"W_Decay_Virtual")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,W_Decay_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  Default_Reader reader;
  if (reader.GetValue<int>("EXTRAXS_W_Decay_Virtual",0) != 1) return NULL;
  // PRINT_VAR(pi.m_fi.m_nloqcdtype << " , " << pi.m_fi.NLOType());
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born) return NULL;
  // if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  // if (pi.m_fi.m_nloewtype&nlo_type::lo) {
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[2].IsLepton() && fl[3].IsLepton() && fl[3]!=fl[2].Bar() && fl[3].LeptonFamily()==fl[2].LeptonFamily() && 
      fl[0].IsQuark()  && fl[1]!=fl[0].Bar() && fl[0].IsUptype()!=fl[1].IsUptype() && fl[0].IsDowntype()!=fl[1].IsDowntype() &&
      (fl[0].IntCharge()+fl[1].IntCharge()==fl[2].IntCharge()+fl[3].IntCharge())) {
    if (pi.m_maxcpl[1]==2 && pi.m_mincpl[1]==2) {
      return new W_Decay_Virtual(pi, fl);
    }
  }
  return NULL;
}


