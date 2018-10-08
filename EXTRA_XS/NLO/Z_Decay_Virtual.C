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
  class Z_Decay_Virtual : public ME2_Base {
    double m_fac;
    int                m_spins[4],m_ew, m_ew_scheme, m_qed, m_lo, m_own;
    Complex m_yqq, m_yll, m_cL_q, m_cL_l, m_cR_q, m_cR_l;
    Flavour    m_flavs[4];
    double m_m2, m_m2p, m_Qf, m_Qfp, m_If, m_Ifp,m_smod;
    double m_e, m_alpha, m_deltas, m_omega, m_mu2;
    double chi1,chi2,qedterm,intterm,Zterm,m_kappa,qf,qe,vf,af,ve,ae,sin2tw,mass;
    double m_quark_masses[5];
    double MW2;
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
    Complex    m_cL[2];
    Complex    m_cR[2];
    METOOLS::DivArrC One,Zero,m_res;
    ATOOLS::Vec4D m_p1,m_p2;
    double m_m_1,m_m_2,m_x1,m_x2,m_xx1,m_xx2;
    std::string widthscheme;
  public:
    Z_Decay_Virtual(const Process_Info& pi, const Flavour_Vector& flavs);

    ~Z_Decay_Virtual();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    Complex InfraredSubtractedME_0_0(const int& b);
    METOOLS::DivArrC InfraredSubtractedME_0_1(const int& b);
    Complex GetBeta_0_0(const int& b1, const int& b2);
    METOOLS::DivArrC GetBeta_0_1(const int& b1, const int& b2);
    METOOLS::DivArrC FAa(const double& s);
    METOOLS::DivArrC FA1(const double& s);
    METOOLS::DivArrC FA3(const double& s);
    METOOLS::DivArrC FV2(const double& s);
    METOOLS::DivArrC FZa(const double& s);
    METOOLS::DivArrC FWa(const double& s);
    METOOLS::DivArrC FWabar(const double& s);
    METOOLS::DivArrC FWn(const double& s);
    METOOLS::DivArrC FWnbar(const double& s);
    METOOLS::DivArrC CT_L(const int& b);
    METOOLS::DivArrC CT_R(const int& b);
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

Z_Decay_Virtual::Z_Decay_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
  ME2_Base(pi, flavs) {
  Data_Reader reader(" ",";","#","=");
  // 0: Only QED plus soft real, 1: weak + QED + soft real, 2: weak (no soft real)
  double ew_prelim = reader.GetValue<int>("YFS_EW_CORRECTIONS",0);
  // EW input scheme: 1: alpha(MZ), 2: alpha(0), 3: Gmu
  m_ew_scheme = reader.GetValue<int>("YFS_EW_SCHEME",2);
  // Turn qed on or off
  m_qed = 1;
  // Check LO if necessary
  m_lo = reader.GetValue<int>("YFS_EXTRAXS_LO",0); 
  if (ew_prelim == 2) m_qed = 0;
  if (ew_prelim > 0) {
    m_ew = 1;
  }
  else m_ew = 0;
  if (m_ew == 0 && m_qed == 0) msg_Out() << "Turned every correction off. This should not happen!\n";
  msg_Debugging() << m_ew << " , " << m_qed << " , " << m_lo << " , " << ew_prelim;
  // Read in width scheme - default is CMS, comparison against paper is in fixed scheme
  widthscheme=reader.GetValue<string>("WIDTH_SCHEME","CMS");
  // Infrared cutoff as used in WZGRAD
  // Conversion: omega_YFS = 0.5*deltas*shat
  m_deltas = reader.GetValue<double>("YFS_DELTAS",0.001);
  Complex I = Complex(0.,1.);
  // Set masses and coupligns for Born
  sin2tw = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  mass     = flavs[2].Mass();
  qe       = flavs[0].Charge();
  qf       = flavs[2].Charge();
  ae       = flavs[0].IsoWeak();      
  af       = flavs[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  m_kappa  = 1./(4.*sin2tw*(1.-sin2tw));


  // Set quark masses to relevant values (useful for test with massles quarks at LO)
  m_quark_masses[0] = (Flavour(kf_d).Mass()==0.?0.06984:Flavour(kf_d).Mass());
  m_quark_masses[1] = (Flavour(kf_u).Mass()==0.?0.06983:Flavour(kf_u).Mass());
  m_quark_masses[2] = (Flavour(kf_s).Mass()==0.?0.15:Flavour(kf_s).Mass());
  m_quark_masses[3] = (Flavour(kf_c).Mass()==0.?1.2:Flavour(kf_c).Mass());
  m_quark_masses[4] = (Flavour(kf_b).Mass()==0.?4.6:Flavour(kf_b).Mass());
  

  // Set boson masses, widths and top mass and width
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
  // if (MODEL::s_model->ScalarNumber("WidthScheme")) {
  muW2 = MW*(MW-I*GW);
  muZ2 = MZ*(MZ-I*GZ);
  muH2 = MH*(MH-I*GH);
  mut2 = mt*(mt-I*Flavour(kf_t).Width());
  // Set mixing angle
  if (widthscheme == "CMS") {
    m_cW2=muW2/muZ2;
    m_sW2=1.-m_cW2;
  }
  else if (widthscheme == "Fixed") {
    // This is the scheme used for comparisons
    m_cW2=MW2/MZ2;
    m_sW2=1.-m_cW2;
  }
  m_sW = sqrt(m_sW2);
  m_cW = sqrt(m_cW2);
  m_cW4 = m_cW2*m_cW2;
  m_sW4 = m_sW2*m_sW2;
  m_e = sqrt(4.*M_PI*AlphaQED());
  m_alpha = AlphaQED();
  // Set flavours for decay MEs
  m_flavs[0]  = Flavour(kf_Z);
  if (!flavs[2].IsAnti()) {
    m_flavs[1] = flavs[2]; 
    m_flavs[2] = flavs[3]; 
  }
  else {
    m_flavs[2] = flavs[2]; 
    m_flavs[1] = flavs[3]; 
  }
  // left/right-handed couplings to be used in X-fucntions
  m_cL[0] = -I*m_e*m_flavs[1].Charge();
  m_cR[0] = -I*m_e*m_flavs[1].Charge();
  m_cL[1] = I*m_e/(2.*m_sW*m_cW)*(2.*m_flavs[1].IsoWeak()
				 -2.*m_flavs[1].Charge()*m_sW2);
  m_cR[1] = I*m_e/(2.*m_sW*m_cW)*(-2.*m_flavs[1].Charge()*m_sW2);

  // Determine isospin partner of lepton
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

  // Set masses, charges for use in EW corrections
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

Z_Decay_Virtual::~Z_Decay_Virtual() {
}

double Z_Decay_Virtual::operator()(const Vec4D_Vector& momenta) {
  // Boost momenta into CMS
  Poincare CMS = Poincare(momenta[0]+momenta[1]);
  Vec4D temp;
  temp = momenta[0]+momenta[1];
  CMS.Boost(temp);
  m_moms[0] = temp;
  temp = momenta[2];
  CMS.Boost(temp);
  m_moms[1] = temp;
  temp = momenta[3];
  CMS.Boost(temp);
  m_moms[2] = temp;
  
  // Calculate Born terms: qedterm ~ |M^gamma|^2, Zterm ~ |M^Z|^2, intterm ~ interference
  double sum = 0.;
  double s(0.),t(0.);
  s=(momenta[0]+momenta[1]).Abs2();
  t=(momenta[0]-momenta[2]).Abs2();
  chi1  = m_kappa * s * (s-MZ2)/(sqr(s-MZ2) + GZ2*MZ2);
  chi2  = sqr(m_kappa * s)/(sqr(s-MZ2) + GZ2*MZ2);
  qedterm = (1+sqr(1.+2.*t/s)) * sqr(qf*qe);
  intterm = (1+sqr(1.+2.*t/s)) * 2.*(qf*qe*vf*ve) * chi1 + (1.+2.*t/s) * 4. * qe*qf*ae*af * chi1;
  Zterm = (1+sqr(1.+2.*t/s)) * (ae*ae+ve*ve) * (af*af+vf*vf) * chi2 + (1.+2.*t/s) * 8. * ae*ve*af*vf * chi2;
  // Check Born if necessary
  if (m_lo == 1) {
    m_res.Finite()=qedterm+intterm+Zterm;
    return sqr(4.*M_PI*m_alpha*CouplingFactor(0,1))*1./3.*m_res.Finite().real();
  }

  // Determine YFS form factor contribution
  DivArrC real = (m_qed?m_alpha*CouplingFactor(0,1)*YFS_Form_Factor()*One:Zero);
  // Reweight to include virtual corrections - note: includes YFS form factor
  DivArrC virtqed = (GetBeta_0_1(0,0)/GetBeta_0_0(0,0)+real)*qedterm;
  DivArrC virtint = ((GetBeta_0_1(0,1)+GetBeta_0_1(1,0))/(GetBeta_0_0(1,0)+GetBeta_0_0(0,1))+real)*intterm;
  DivArrC virtZ = (GetBeta_0_1(1,1)/GetBeta_0_0(1,1)+real)*Zterm;

  // Compare to own implementation.
  m_own = 1;

  if (m_qed == 1 && m_ew == 0 && m_lo == 0 && m_own == 0) {
    // finite term as in WZGRAD - sanity check only
    m_res.Finite()=m_alpha*CouplingFactor(0,1)/M_PI*(2.*(log(s/m_m2)-1.)*log(m_deltas)+3./2.*log(s/m_m2)-2.-sqr(M_PI)/6.)*(qedterm+intterm+Zterm);
  }
  else {
    // UV, IR, IR2 terms not strictly needed but nice for checking divergence cancellation ... 
    // 1/epsUV
    m_res.IR()=(virtqed+virtint+virtZ).UV().real();
    // 1/epsIR
    m_res.IR()=(virtqed+virtint+virtZ).IR().real();
    // 1/epsIR2
    m_res.IR2()=(virtqed+virtint+virtZ).IR2().real();
    // finite
    m_res.Finite()=(virtqed+virtint+virtZ).Finite().real();
  }
  return sqr(4.*M_PI*m_alpha*CouplingFactor(0,1))*1./3.*m_res.Finite().real();
}


Complex Z_Decay_Virtual::InfraredSubtractedME_0_0(const int& b) {
  // b = 0: photon exchanged, b = 1: Z exchanged
  // get M_0^0 for either gamma^* or Z in initial state
  // eps is the same as gamma^* is vector of mass q^2
  Vec4C epsZ = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  return XYZ.X(1,m_spins[1],epsZ,2,m_spins[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
}

DivArrC Z_Decay_Virtual::InfraredSubtractedME_0_1(const int& b) {
  // Set up ingredients for ME calculation
  Vec4C epsV = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  double s((m_moms[1]+m_moms[2]).Abs2());
  double m2(sqr(m_flavs[1].Mass()));
  double mu2(s);
  m_mu2 = s;
  // Set up EW corrections as corrections to cR,cL
  // If only QED corrections, all weak form factors and QED contributions in CTs are neglected
  DivArrC term(0.,0.,0.,0.,0.,0.);
  if (b == 0) {
    // photon exchange
    // Collect form factors
    DivArrC cRDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(-Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))
				     *(
				       (m_qed?m_Qf/4.*FAa(s):Zero)
				       +(m_ew?m_Qf/(4.*m_sW2*m_cW2)*pow(-m_sW2*m_Qf,2.)*FZa(s):Zero)
				       )
				     +CT_R(b)
				     );    
    DivArrC cLDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(-Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))
				     *(
				       (m_qed?m_Qf/4.*FAa(s):Zero)
				       +(m_ew?(m_Qf/(4.*m_sW2*m_cW2)*pow(m_If-m_sW2*m_Qf,2.)*FZa(s)
					       +m_Qfp/(8.*m_sW2)*FWa(s)
					       -m_If/(4.*m_sW2)*FWn(s)):DivArrC(0.,0.,0.,0.,0.,0.))
				       )
				     +CT_L(b)
				     );  
    // Contributions proportional to Y(1,2,cR,cL) 
    DivArrC cRDivArrY = (m_qed?m_alpha*CouplingFactor(0,1)/M_PI*(m_e*sqrt(CouplingFactor(0,1))
								 *m_Qf/4.*FV2(s)*(epsV*(m_moms[1]-m_moms[2]))):Zero);
    DivArrC cLDivArrY = (m_qed?m_alpha*CouplingFactor(0,1)/M_PI*(m_e*sqrt(CouplingFactor(0,1))
				      *m_Qf/4.*FV2(s)*(epsV*(m_moms[1]-m_moms[2]))
								 ):Zero);

    // Calculate IR factor B - not needed if only weak parts are checked
    DivArrC B(0.,0.,0.,0.,0.,0.);
    // m2*C_0(m2,m2,0,0,m2,m2) contains an IR divergence necessary to cancel all divergences even 
    // when m2 = 0 (C_0 ~ 1/m2)
    if (m_qed) {
      if (m2 != 0) {
	B +=  m_alpha*CouplingFactor(0,1)/M_PI*(0.5*(s-2.*m2)*C_0(m2,m2,s,m2,0.,m2,mu2)
						+m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
						+0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2)))
	  *XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
      }
      else if (m2 == 0) {
	B +=  m_alpha*CouplingFactor(0,1)/M_PI*(0.5*(s-2.*m2)*C_0(m2,m2,s,0.,m2,m2,mu2)
						+0.5*DivArrC(0.,1.,0.,0.,0.,0.)
						+0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2)))
	  *XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
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
  else if (b == 1) {
    // Z exchange
    // Assemble form factors
    DivArrC cRDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW)
				     *(
				       (m_qed?(-m_sW2*m_Qf)/4.*FAa(s)
				     	-m_If/8.*FA1(s):Zero)
				       +(m_ew?1./(4.*m_sW2*m_cW2)*pow(-m_sW2*m_Qf,3.)*FZa(s):Zero)
				       )
				     +
							 CT_R(b)
				     );    
    DivArrC cLDivArr = m_alpha*CouplingFactor(0,1)/M_PI*(Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW)
				     *(
				       (m_qed?(m_If-m_sW2*m_Qf)/4.*FAa(s)
				     	+m_If/8.*FA1(s):Zero)
				       +(m_ew?(1./(4.*m_sW2*m_cW2)*pow(m_If-m_sW2*m_Qf,3.)*FZa(s)
				     	       +1./(8.*m_sW2)*((m_Ifp-m_sW2*m_Qfp)*FWa(s)+m_Ifp/2.*FWabar(s))
				     	       -m_If*m_cW2/(4.*m_sW2)*(FWn(s)+FWnbar(s))):Zero)
				       )
				     +
							 CT_L(b)
				     );  
    // Terms proportional to Y(1,2,cR,cL)
    DivArrC cRDivArrY = (m_qed?m_alpha*CouplingFactor(0,1)/M_PI*(m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW)
								 *((m_If-2.*m_sW2*m_Qf)/8.*FV2(s)*(epsV*(m_moms[1]-m_moms[2]))
								   -m_If/8.*FA3(s)*(epsV*(m_moms[1]+m_moms[2])))):Zero);
    DivArrC cLDivArrY = (m_qed?m_alpha*CouplingFactor(0,1)/M_PI*(m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW)
								 *((m_If-2.*m_sW2*m_Qf)/8.*FV2(s)*(epsV*(m_moms[1]-m_moms[2]))
								   +m_If/8.*FA3(s)*(epsV*(m_moms[1]+m_moms[2])))):Zero);
    // Calculate IR factor B - not needed if only weak parts are checked
    DivArrC B(0.,0.,0.,0.,0.,0.);
    // m2*C_0(m2,m2,0,0,m2,m2) contains an IR divergence necessary to cancel all divergences even 
    // when m2 = 0 (C_0 ~ 1/m2)
    if (m_qed) {
      if (m2 != 0) {
	B +=  m_alpha*CouplingFactor(0,1)/M_PI*(0.5*(s-2.*m2)*C_0(m2,m2,s,m2,0.,m2,mu2)
						+m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
						+0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2)))
	  *XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
      }
      else if (m2 == 0) {
	B +=  m_alpha*CouplingFactor(0,1)/M_PI*(0.5*(s-2.*m2)*C_0(m2,m2,s,0.,m2,m2,mu2)
						+0.5*DivArrC(0.,1.,0.,0.,0.,0.)
						+0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2)))
	  *XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
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
}

Complex Z_Decay_Virtual::GetBeta_0_0(const int& b1, const int& b2) {
  // Assemble Born decay ME squared
  // b1: boson exchanged in amplitude, b2: boson exchange in conjugated amplitude
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
        m_spins[0] = k;
        m_spins[1] = j;
        m_spins[2] = i;
        Complex M_0_0_b1 = InfraredSubtractedME_0_0(b1);
        Complex M_0_0_b2 = InfraredSubtractedME_0_0(b2);
        sum = sum + (M_0_0_b1*conj(M_0_0_b2))// .real()
	  ;
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  return sum;
}

DivArrC Z_Decay_Virtual::GetBeta_0_1(const int& b1, const int& b2) {
  // Assemble virtual decay ME squared
  // b1: boson exchanged in tree amplitude, b2: boson exchange in virtual amplitude
  DivArrC sum(0.,0.,0.,0.,0.,0.);
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
	m_spins[0] = k;
	m_spins[1] = j;
	m_spins[2] = i;
	Complex M_0_0 = InfraredSubtractedME_0_0(b1);
	DivArrC M_0_1 = InfraredSubtractedME_0_1(b2);
	sum = sum + (M_0_0*conj(M_0_1)+conj(M_0_0)*M_0_1);
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  return sum;
}



// Form factors for virtual amplitude (numbers refer to equations in
// Bardin:1999ak)
DivArrC Z_Decay_Virtual::B_ff(const double& p2, const Complex& m12,
					 const Complex& m22) {
  // Helper function
  return B_0(p2,m12,m22,p2) - B_0(real(m12),m12,0.,p2);
}

// Photon exchange (5.568)
DivArrC Z_Decay_Virtual::FAa(const double& s) {
  return -2.*(s-2*m_m2)*C_0(m_m2,m_m2,s,m_m2,0.,m_m2,m_mu2) + B_0(s,m_m2,m_m2,m_mu2)
    - 4.*B_ff(s,m_m2,m_m2) - 2.*One;
}

// Additional mass dependent photon pieces (5.567), (5.588)
DivArrC Z_Decay_Virtual::FA1(const double& s) {
  return -2.*s*m_m2/(s*(4.*m_m2-s)/4.)*B_ff(s,m_m2,0.);
}

DivArrC Z_Decay_Virtual::FA3(const double& s) {
  return -sqrt(m_m2)/(s*(4.*m_m2-s)/4.)*((4.*m_m2-3.*s)/2.*B_ff(s,m_m2,0.)
					     + (4.*m_m2 - s)*One);
}

DivArrC Z_Decay_Virtual::FV2(const double& s) {
  return s*sqrt(m_m2)/(s*(4.*m_m2-s)/2.)*B_ff(s,m_m2,m_m2);
}

// Z exchange (5.571)
DivArrC Z_Decay_Virtual::FZa(const double& s) {
  Complex z(-s/muZ2);
  return 
    2.*muZ2/z*pow((1.-z),2.)*C_0(0.,0.,s,0.,muZ2,0.,m_mu2)
    + B_0(s,0.,0.,m_mu2)
    + (2./z-4.)*(B_0(s,0.,0.,m_mu2)-B_0(0.,0.,muZ2,m_mu2))
    -2.*One;
}

// W exchange (5.577)
DivArrC Z_Decay_Virtual::FWa(const double& s) {
  Complex w(-s/muW2);
  Complex wf = m_m2/muW2;
  Complex zf = m_m2/muZ2;
  Complex wfp = m_m2p/muW2;
  Complex zfp = m_m2p/muZ2;
  Complex beta2 = 1. - wfp;
  Complex kappa = -0.5*beta2*(3.-beta2)*muW2/s; 
  return -(-2.*beta2*kappa + 3. + pow(beta2,2.) -2.*w)*muW2*C_0(0.,0.,s,m_m2p,muW2,m_m2p,m_mu2)
    +(3.-beta2)/2.*B_0(s,m_m2p,m_m2p,m_mu2)
    +2.*(kappa-2.)*(B_0(s,m_m2p,m_m2p,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
    -(2.+0.5*wfp)*One;
}

// (5.593)
DivArrC Z_Decay_Virtual::FWabar(const double& s) {
  Complex w(-s/muW2);
  Complex wf = m_m2/muW2;
  Complex zf = m_m2/muZ2;
  Complex wfp = m_m2p/muW2;
  Complex zfp = m_m2p/muZ2;
  Complex beta2 = 1. - wfp;
  Complex kappa = -0.5*beta2*(3.-beta2)*muW2/s; 
  return wfp*(-(pow(beta2,2.)/w-2.)*muW2*C_0(0.,0.,s,m_m2p,muW2,m_m2p,m_mu2)
	      - 0.5*B_0(s,m_m2p,m_m2p,m_mu2)
	      - beta2/w*(B_0(s,m_m2p,m_m2p,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
	      + 0.5*One);
}

// Triple W vertex (5.581)
DivArrC Z_Decay_Virtual::FWn(const double& s) {
  Complex w(-s/muW2);
  Complex wf = m_m2/muW2;
  Complex zf = m_m2/muZ2;
  Complex wfp = m_m2p/muW2;
  Complex zfp = m_m2p/muZ2;
  Complex beta2 = 1. - wfp;
  Complex kappa = -0.5*beta2*(3.-beta2)*muW2/s; 
  return -(-2.*beta2*kappa+3.+pow(beta2,2.))*muW2*C_0(0.,0.,s,muW2,m_m2p,muW2,m_mu2)
    -2.*(kappa-2.)*(B_0(s,muW2,muW2,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
    -(3.+0.5*wfp)*B_0(s,muW2,muW2,m_mu2)
    -0.5*wfp*One;
}

// (5.596)
DivArrC Z_Decay_Virtual::FWnbar(const double& s) {
  Complex w(-s/muW2);
  Complex wf = m_m2/muW2;
  Complex zf = m_m2/muZ2;
  Complex wfp = m_m2p/muW2;
  Complex zfp = m_m2p/muZ2;
  Complex beta2 = 1. - wfp;
  Complex kappa = -0.5*beta2*(3.-beta2)*muW2/s; 
  return 0.5/m_cW2*wfp*(-(pow(beta2,2.)/w+4.-wfp)*muW2*C_0(0.,0.,s,muW2,m_m2p,muW2,m_mu2)
			+0.5*B_0(s,muW2,muW2,m_mu2)
			+beta2/w*(B_0(s,muW2,muW2,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
			+0.5*One);
}


DivArrC Z_Decay_Virtual::CT_R(const int& b)
{
  // Counterterm contribution to right-handed coupling (proportional to PR) 
  // b = 0: eey counterterm, b = 1: eeZ counterterm
  if (b == 0) {
    Complex pref = Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1));
    Complex gR = -m_sW/m_cW*m_Qf; // Note: different from the below gR due to prefactor ...
    if (!m_ew) return -m_Qf*pref*(1./2.*dZAA()+real(dZe())+dZfermR(m_m2,m_m2p,m_Qf,m_If));
    else return -m_Qf*pref*(
			    1./2.*dZAA()
			    +dZe()
			    +dZfermR(m_m2,m_m2p,m_Qf,m_If))
	   +1./2.*pref*gR*dZZA();
    
  }
  else if (b == 1) {
    Complex pref = Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW);
    Complex gR = -m_sW2*m_Qf;
    if (!m_ew) return pref*gR*(dZfermR(m_m2,m_m2p,m_Qf,m_If));
    else return pref*gR*(
			 1./2.*dZZZ()
			 +dZe()
			 -dcw()/m_sW2
			 +dZfermR(m_m2,m_m2p,m_Qf,m_If))
	   -1./2.*Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))*m_Qf*dZAZ();
  }
}

DivArrC Z_Decay_Virtual::CT_L(const int& b)
{
  // Counterterm contribution to left-handed coupling (proportional to PL)
  // b = 0: eey counterterm, b = 1: eeZ counterterm
  if (b == 0) {
    Complex pref = Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1));
    Complex gL = (m_If-m_sW2*m_Qf)/(m_sW*m_cW);
    if (!m_ew) return -m_Qf*pref*(dZfermL(m_m2,m_m2p,m_Qf,m_If));
    else return -m_Qf*pref*(
			   1./2.*dZAA()
			 +dZe()
			 +dZfermL(m_m2,m_m2p,m_Qf,m_If))
	   +1./2.*pref*gL*dZZA();
  }
  else if (b == 1) {
    Complex pref = Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(m_sW*m_cW);
    Complex gL = m_If-m_sW2*m_Qf;
    if (!m_ew) return pref*gL*(dZfermL(m_m2,m_m2p,m_Qf,m_If));
    else return pref*gL*(
			 1./2.*dZZZ()
			 +dZe()
			 +dZfermL(m_m2,m_m2p,m_Qf,m_If))
	   +pref*dcw()*(m_If*(m_cW2-m_sW2)/m_sW2+m_Qf)
	   -1./2.*Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))*m_Qf*dZAZ();
  }
}






// Calculate YFS form factor - copied and adapted from PHOTONS
double Z_Decay_Virtual::YFS_Form_Factor() {
  double s = (m_moms[1]+m_moms[2]).Abs2();
  double mu2 = s;
  // Set IR cutoff
  m_omega = m_deltas*sqrt(s)/2.;
  // Set dipole momenta and masses, and order by energy (not really necessary here)
  m_p1 = m_moms[1];
  m_p2 = m_moms[2];
  m_m_1 = m_flavs[1].Mass();//m_moms[1].Abs(); // Do not use Abs() here as this can create numerical problems with electrons
  m_m_2 = m_flavs[2].Mass();//m_moms[2].Abs();
  if (m_moms[2][0] >= m_moms[1][0]) {
    std::swap(m_p1,m_p2);
    std::swap(m_m_1,m_m_2);
  }
  // Calculate and store the x_i and x'_i used in form factor calculation
  m_x1  = - (m_p1.Abs2() - m_p2.Abs2()
	     + 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
    / (m_p1-m_p2).Abs2();
  m_x2  = - (m_p1.Abs2() - m_p2.Abs2()
	     - 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
    / (m_p1-m_p2).Abs2();
  if (abs(1.+m_x1) < 1E-10) {
    msg_Error()<<METHOD<<"() error: case should not appear !!!\n";
  }
  if (abs(1.+m_x2) < 1E-10) {
    // expansion for (p1-p2)<~0 => x2>~-1
    double p22(m_p2.Abs2());
    double x(m_p1.Abs2()*m_p2.Abs2()/sqr(m_p1*m_p2));
    double r(m_p1*m_p2);
    m_x2  = - 1. + (2.*p22-r*(-0.5*x+0.125*x*x-0.0625*x*x*x))
      / (m_p1-m_p2).Abs2();
  }
  m_xx1 = - (m_p1.Abs2() - m_p2.Abs2()
	     + 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
    / (m_p1+m_p2).Abs2();
  m_xx2 = - (m_p1.Abs2() - m_p2.Abs2()
	     - 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
    / (m_p1+m_p2).Abs2();
  // Assemble form factor. Added pi^2 in to match WZGRAD
  return 1./M_PI*(log(m_p1[0]*m_p2[0]/sqr(m_omega)) + 0.5*(m_p1*m_p2)*IntP1()
		  - 0.5*(m_p1*m_p2)*IntE() + 0.25*IntP2() + G(1.) + G(-1.)
		  - (m_p1*m_p2)*IntG() + sqr(M_PI));
}

double Z_Decay_Virtual::IntP1() {
  // first p^2 integral (contains 1/lambda)
  double A(0.);
  if (m_xx1*m_xx2 >= 0.) 
    A = 8.*(M_PI*M_PI)/((m_xx2-m_xx1)*(m_p1+m_p2).Abs2());
  double B = 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
    *(log(abs(m_x1))*(DiLog((m_x1-1.)/m_x1)-DiLog((m_x1+1.)/m_x1))
      -log(abs(m_x2))*(DiLog((m_x2-1.)/m_x2)-DiLog((m_x2+1.)/m_x2)));
  return (A+B);
}

double Z_Decay_Virtual::IntE() {
  // E1 = E2
  // integral involving energy
  if (abs(m_p1[0]-m_p2[0]) < 1E-6) {
    return 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
           *log((m_p1[0]+m_p2[0])/(2*m_omega))
           *log(abs(((1.-m_x1)*(1.+m_x2))/((1.+m_x1)*(1.-m_x2))));
  }
  else {
    msg_Out() << METHOD << "Leptons don't have the same energy?\n";
    return 0.;
  }
}

// int dx ln (px'²/m1m2)                                                  (C.55)
double Z_Decay_Virtual::IntP2() {
  return 2.*log((m_p1+m_p2).Abs2()/(4.*m_m_1*m_m_2))
    + log(abs((1.-m_xx1*m_xx1)*(1.-m_xx2*m_xx2)))
    - m_xx1*log(abs((1.-m_xx1)/(1.+m_xx1)))
    - m_xx2*log(abs((1.-m_xx2)/(1.+m_xx2))) - 4.;
}

// G(x)
double Z_Decay_Virtual::G(double x) {
  Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
  double b  = CalculateBeta(px);
  double r(0.);
  if (b == 0.)      r = 1.-LOG_2;
  else if (b == 1.) r = 0.;
  else              r = ((1.-b)/(2.*b)*log((1.+b)/(1.-b))+log((1.+b)/2.));
  return r;
}

// int dx G(x)/px²                                                        (C.88)
double Z_Decay_Virtual::IntG() {
  // if dipole in its CMS
  if ((Vec3D(m_p1)+Vec3D(m_p2)).Abs() < 1E-3) {
    // same mass or both nearly massless or both of nearly same beta
    if ((abs(m_m_1-m_m_2) < 1E-6) || 
        ((1.-CalculateBeta(m_p1) < 5E-3) && (1.-CalculateBeta(m_p2) < 5E-3)) ||
        ((CalculateBeta(m_p1)-CalculateBeta(m_p2))
          /(CalculateBeta(m_p1)+CalculateBeta(m_p2)) < 5E-3)) {
      double E = m_p1[0];
      double b = CalculateBeta(m_p1);
      double r = 1./(b*E*E)*(1./2.*sqr(log((1.+b)/2.)) + LOG_2*log(1.+b)
                             - 1./2.*sqr(LOG_2) - 1./2.*sqr(log(1.+b))
                             + DiLog((1.-b)/2.) - DiLog((1.+b)/2.)
                             + DiLog(b) - DiLog(-b));
      return r;
    }
  }
  else {
    msg_Out() << METHOD << "Leptons not in dipole CMS or of different mass?\n";
  }
}

// Helper function
double Z_Decay_Virtual::CalculateBeta(const Vec4D& p) {
  return (Vec3D(p).Abs()/p[0]);
}








// Functions to determine EW counterterms - taken from Denner:1991kt
DivArrC Z_Decay_Virtual::dZfermL(const double& m2, const double& m2p,
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

DivArrC Z_Decay_Virtual::dZfermR(const double& m2, const double& m2p,
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

DivArrC Z_Decay_Virtual::dm(const double& m2, const double& m2p,
			    const double& Qf, const double& If) {
  // Fermion mass counterterm
  return sqrt(m2)/2.*(Sigma_ferm_L(m2,m2,m2p,Qf,If) 
		      + Sigma_ferm_R(m2,m2,m2p,Qf,If) 
		      + 2.*Sigma_ferm_S(m2,m2,m2p,Qf,If));
}

// Fermion self energies
DivArrC Z_Decay_Virtual::Sigma_ferm_R(const double& p2,
				      const double& m2, const double& m2p,
				      const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf);
  return (m_qed?-1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One):Zero)
    +(m_ew?(-1./4.*pow(gplus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	    -1./(16.*m_sW2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2/muW2*B_1(p2,m2p,muW2,m_mu2)):Zero);
}

DivArrC Z_Decay_Virtual::Sigma_ferm_L(const double& p2,
				      const double& m2, const double& m2p,
				      const double& Qf, const double& If)
{
  Complex gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One):Zero)
    +(m_ew?(-1./4.*pow(gminus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	    -1./(16.*m_sW2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*((2.+m2p/muW2)*B_1(p2,m2p,muW2,m_mu2)+1.*One)):Zero);
}


DivArrC Z_Decay_Virtual::Sigma_ferm_S(const double& p2,
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
DivArrC Z_Decay_Virtual::dSigma_ferm_R(const double& p2,
				       const double& m2, const double& m2p,
				       const double& Qf, const double& If)
{
  Complex gplus(-m_sW/m_cW*Qf);
  return (m_qed?-1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2):Zero)
    +(m_ew?(-1./4.*pow(gplus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sW2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*m2/muW2*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}

DivArrC Z_Decay_Virtual::dSigma_ferm_L(const double& p2,
				       const double& m2, const double& m2p,
				       const double& Qf, const double& If)
{
  Complex gminus((If-m_sW2*Qf)/(m_sW*m_cW));
  return (m_qed?-1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2):Zero)
    +(m_ew?(-1./4.*pow(gminus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sW2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sW2)*(2.+m2p/muW2)*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}


DivArrC Z_Decay_Virtual::dSigma_ferm_S(const double& p2,
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
DivArrC Z_Decay_Virtual::Sigma_W(const double& p2)
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
    double m2d(sqr(m_quark_masses[2*i])), m2u(((i==2)?mt2:sqr(m_quark_masses[2*i+1])));
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


DivArrC Z_Decay_Virtual::dSigma_W(const double& p2)
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
    double m2d(sqr(m_quark_masses[2*i])), m2u(((i==2)?mt2:sqr(m_quark_masses[2*i+1])));
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
DivArrC Z_Decay_Virtual::Sigma_ZZ(const double& p2)
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
DivArrC Z_Decay_Virtual::dSigma_ZZ(const double& p2)
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
DivArrC Z_Decay_Virtual::Sigma_H(const double& p2)
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
DivArrC Z_Decay_Virtual::dSigma_H(const double& p2)
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
DivArrC Z_Decay_Virtual::dSigma_GamGam(const double& p2)
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
DivArrC Z_Decay_Virtual::dSigma_GamGam0()
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
     +1./3.*(2.*mb2*B_0(0.,mb2,mb2,m_mu2) - (2.*mb2+muZ2)*B_0(MZ2,mb2,mb2,m_mu2))) 
     +muZ2*(1.*One+11./27.*One);
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
DivArrC Z_Decay_Virtual::Sigma_ZA(const double& p2)
{
  if (m_ew == 0) return Zero;
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:sqr(m_quark_masses[i]))),ch(flav.Charge());
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
DivArrC Z_Decay_Virtual::dZW()
{
  // return -real(dSigma_W(muW2));
  return -dSigma_W(MW2);
}

DivArrC Z_Decay_Virtual::dZZZ()
{
  // return -real(dSigma_ZZ(MZ2));
  return -dSigma_ZZ(MZ2);
}

DivArrC Z_Decay_Virtual::dZAA()
{
  if (m_ew_scheme == 1 || m_ew_scheme == 3) return -dSigma_GamGam(MZ2);
  else if (m_ew_scheme == 2) return -dSigma_GamGam0();
}

DivArrC Z_Decay_Virtual::dZAZ()
{
  if (widthscheme == "Fixed") return -Complex(1.,0.)*2.*real(Sigma_ZA(MZ2))/MZ2;
  else if (widthscheme == "CMS") return -2.*Sigma_ZA(MZ2)/MZ2 + (muZ2/MZ2 - 1.)*dZZA();
}

DivArrC Z_Decay_Virtual::dZZA()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*2.*Sigma_ZA(0.)/muZ2;
  else if (widthscheme == "CMS") return 2.*Sigma_ZA(0.)/muZ2;
}

DivArrC Z_Decay_Virtual::dZH()
{
  return -dSigma_H(MH2);
}

// Mass Counterterms
DivArrC Z_Decay_Virtual::dMW2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_W(MW2));
  else if (widthscheme == "CMS") return Sigma_W(MW2) + (muW2 - MW2)*dSigma_W(MW2);
}

DivArrC Z_Decay_Virtual::dMZ2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_ZZ(MZ2));
  else if (widthscheme == "CMS") return Sigma_ZZ(MZ2) + (muZ2-MZ2)*dSigma_ZZ(MZ2);
}

DivArrC Z_Decay_Virtual::dMH2()
{
  if (widthscheme == "Fixed") return Complex(1.,0.)*real(Sigma_H(MH2));
  else if (widthscheme == "CMS") return Sigma_H(MH2) + (muH2-MH2)*dSigma_H(MH2);
}

// Counterterm for cos(thetaW)
DivArrC Z_Decay_Virtual::dcw()
{
  if (widthscheme == "Fixed") return 1./2.*(dMW2()/MW2-dMZ2()/MZ2);
  else if (widthscheme == "CMS") return 1./2.*(dMW2()/muW2-dMZ2()/muZ2);
}

// Charge renormalization counterterm
DivArrC Z_Decay_Virtual::dZe()
{
  switch (m_ew_scheme) {
  case 1: {
    // alpha(MZ)
    if (widthscheme == "Fixed") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/MZ2;  
    if (widthscheme == "CMS") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/muZ2;  
    break;
  }
  case 2: {
    // alpha(0)
    return 1./2.*dSigma_GamGam0() - m_sW/m_cW*Sigma_ZA(0.)/muZ2;
    break;
  }
  case 3: {
    // Gmu 
    return -m_cW2/m_sW2*dcw() + 1./2.*(Sigma_W(MW2)-Sigma_W(0.))/muW2 
      - 1./(m_cW*m_sW)*Sigma_ZA(0.)/muZ2 - One*1./8.*1./m_sW2*(6.+(7.-4.*m_sW2)/(2.*m_sW2)*log(m_cW2));
    break;
  }
  default: { 
    if (widthscheme == "Fixed") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/MZ2; 
    if (widthscheme == "CMS") return 1./2.*dSigma_GamGam(MZ2) - m_sW/m_cW*Sigma_ZA(0.)/muZ2; 
    break;
  }
  }
}






DECLARE_TREEME2_GETTER(Z_Decay_Virtual,"Z_Decay_Virtual")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,Z_Decay_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  Default_Reader reader;
  // Only use this ME if explicitly asked for this particular ME
  // Prevents DY ME not being used
  if (reader.GetValue<int>("EXTRAXS_Z_Decay_Virtual",0) != 1) return NULL;
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if ((fl[2].IsLepton() && fl[3]==fl[2].Bar() &&
       fl[0].IsQuark()  && fl[1]==fl[0].Bar()) ||   
      (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
       fl[2].IsQuark()  && fl[3]==fl[2].Bar())) {
    if (pi.m_maxcpl[1]==2 && pi.m_mincpl[1]==2) {
      return new Z_Decay_Virtual(pi, fl);
    }
  }
  return NULL;
}


