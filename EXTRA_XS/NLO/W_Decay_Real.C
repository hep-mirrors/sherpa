#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Math/Tensor_Build.H"
#include "ATOOLS/Math/Tensor_Contractions.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include <map>
#include <string>


using namespace ATOOLS;
using namespace METOOLS;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS {

  class W_Decay_Real : public ME2_Base {
  private:
    int                m_spins[4];
    int                m_spins0[3];
    int                m_quarkspins[2];
    int                m_1to2spins[3];
    int                m_1to3spins[4];
    int number_pol,pol_dep, m_counter;
    Flavour    m_flavs[9];
    Flavour    m_quark_flavs[2], m_lepton;
    Flavour    m_born_flavs[3], m_real_flavs[6];
    double m_e, m_alpha, m_deltas, m_omega, m_mu2, chi, born, m_below;
    double chi1,chi2,qedterm,intterm,Zterm,m_kappa,qf,qe,vf,af,ve,ae,sin2tw,mass;
    Complex             m_sW, m_sW2, m_sW4;
    Complex             m_cW, m_cW2, m_cW4;
    double m_theta, m_phi, m_c, m_max, m_u;
    double m_quark_masses[5];
    double MW2, GW2, MW;
    double MZ2, GZ2;
    double m_angle_weight;
    Complex muW2;
    Complex muZ2, m_total_weight, m_total_weight2;
    ATOOLS::Vec4D m_dir;
    ATOOLS::Vec4D m_moms[4];
    ATOOLS::Vec4D m_moms0[3];
    ATOOLS::Vec4D m_quarkmoms[2];
    ATOOLS::Vec4D m_1to2moms[3],m_1to3moms[4];
    double m_masses[3];
    Complex    m_cL;
    Complex    m_cR;
    std::map<std::string, ATOOLS::Histogram * >m_histos;
  public:

    W_Decay_Real(const Process_Info& pi,const Flavour_Vector& fl);
    ~W_Decay_Real();
    Complex InfraredSubtractedME_0_0();
    Complex InfraredSubtractedME_0_0(const int& Wpol);
    Complex ProdME(const int& Wpol);
    Complex RealProdME(const int& Wpol);
    Complex LO_DecayME(const int& Wpol);
    Complex ProdMEconj(const int& Wpol);
    Complex LO_DecayMEconj(const int& Wpol);
    Complex InfraredSubtractedME_1_05();
    Complex InfraredSubtractedME_1_05(const int& Wpol, const int& Gpol);
    Complex GetBeta_0_0();
    Complex GetBeta_1_1();
    double Smod();
    double Smod_Vecs(const ATOOLS::Vec4D& pW, const ATOOLS::Vec4D& pl,const ATOOLS::Vec4D& pnu,const ATOOLS::Vec4D& pg);
    double Smod_Uncorr();
    ATOOLS::Vec4D Generate_One_Photon(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    double Generate_Energy(const double& omega, const double& max);
    void Generate_Photon_Angle(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    void GenerateDipoleAngle(const double& b1, const double& b2);
    void GenerateDipoleAngle(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    void GenerateNullVector();
    double Weight_Jac_M();
    double Weight_Jac_L();
    void DetermineU(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    double DetermineMaxEnergy();
    void Setup_Corrected_Momenta(const ATOOLS::Vec4D_Vector& mom);
    void InitAnalysis();
    void FinishAnalysis(const ATOOLS::Flavour& fl_lep);
    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

W_Decay_Real::W_Decay_Real
(const Process_Info& pi, const Flavour_Vector& fl) : ME2_Base(pi, fl)
{
  Data_Reader reader(" ",";","#","=");
  // deltas: translates to min. photon energy in partonic CMS frame E_min = 0.5*deltas*shat
  m_deltas = reader.GetValue<double>("YFS_DELTAS",0.001);
  string widthscheme=reader.GetValue<string>("WIDTH_SCHEME","CMS");
  // Need to generate lower energies than cutoff, since generation in different frame. 10^(-3) or 10^(-4) should generally be enough
  m_below = reader.GetValue<double>("YFS_BELOW_DELTAS",0.0001);
  // Number of W polarizations - 1. Only used for silly checks.
  number_pol = reader.GetValue<int>("YFS_NUMBER_POL",2);
  // Calculation based on assembly of X- and Y-function (pol_dep=1) or on reweighting of Born with decay ME (pol_dep=0)
  pol_dep = reader.GetValue<int>("YFS_POL_DEPENDENT",0);
  Complex I = Complex(0.,1.);
  sin2tw = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));

  // Set up flavours for born and real amplitude.
  m_born_flavs[0] = Flavour(kf_Wplus,0);
  m_born_flavs[1] = fl[3];
  m_born_flavs[2] = fl[2];
  for (int i=0; i < 3; i++) {
    m_real_flavs[i] = m_born_flavs[i];
  }
  m_real_flavs[3] = Flavour(kf_photon);

  // Helper for analyses
  m_lepton = m_born_flavs[1];
  // Electroweak parameters
  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  MW2 = pow(MW,2.);
  MZ2 = pow(MZ,2.);
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  GZ2 = GZ*GZ;
  // This is the zeroth order decay width
  GW2 = GW*GW;
  muW2 = MW*(MW-I*GW);
  muZ2 = MZ*(MZ-I*GZ);
  if (widthscheme == "CMS") {
    m_cW2=muW2/muZ2;
    m_sW2=1.-m_cW2;
  }
  else if (widthscheme == "Fixed") {
    m_cW2=MW2/MZ2;
    m_sW2=1.-m_cW2;
  }
  m_sW = sqrt(m_sW2);
  m_cW = sqrt(m_cW2);
  m_e = sqrt(4.*M_PI*AlphaQED());
  m_alpha = AlphaQED();

  m_flavs[0] = m_born_flavs[0];
  m_flavs[1] = m_born_flavs[1];
  m_flavs[2] = m_born_flavs[2];
  m_flavs[3] = m_born_flavs[3];
  m_cL = Complex(1.,0.);
  m_cR = Complex(0.,0.);
  InitAnalysis();
  msg->SetPrecision(16); 
  // Helpers for total XS and deviation
  m_total_weight = 0.;
  m_total_weight2 = 0.;
  m_counter = 0.;
}

W_Decay_Real::~W_Decay_Real()
{
  FinishAnalysis(m_lepton);
  Complex XS = m_total_weight/((double)m_counter);
  PRINT_VAR("XS = " << XS);
  PRINT_VAR("sqrt(total_weight2/N) = " << sqrt(m_total_weight2/sqr(((double)m_counter))));
  PRINT_VAR("sqrt((XS^2-total_weight2)/N) = " << sqrt((sqr(XS)-m_total_weight2)/sqr(((double)m_counter))));
}

void W_Decay_Real::InitAnalysis() {
  // Initialize histograms
  m_histos[string("E_lep")]	 = new Histogram(0,0.,55.,55);
  m_histos[string("E_neu")]	 = new Histogram(0,0.,55.,55);
  m_histos[string("Mlgamma")]	 = new Histogram(0,0.,100.,100);
  m_histos[string("E_gamma")]	 = new Histogram(0,0.,55.,55);
  m_histos[string("Mlgamma_small")]	 = new Histogram(0,0.,10.,100);
  m_histos[string("E_gamma_small")]	 = new Histogram(0,0.,10.,100);
  m_histos[string("costheta")]	 = new Histogram(0,-1.,1.,200);
}

void W_Decay_Real::FinishAnalysis(const Flavour& fl_lep) {
  // Finalize histograms
  Histogram * histo;
  string name;
  string lepname = fl_lep.ShellName();
  for (map<string,Histogram *>::iterator 
     hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    histo->MPISync();
    name  = string("Debug_Analysis_")+lepname+string("/")+hit->first+string(".dat");
    PRINT_VAR(name);
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

double W_Decay_Real::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  // New setup for decay only. Use this class only as a vehicle for
  // own calculation Only concerned with m_lgamma and E_i so start off
  // with Born configuration with lepton and neutrino back-to-back
  m_counter++;
  m_1to2moms[0] = Vec4D(MW,0.,0.,0.);

  // Energy cutoff in generation. Take very low value to be most inclusive.
  m_omega = 5.*pow(10.,-4.);

  m_masses[0] = m_1to2moms[0].Abs();
  for (int i(1); i < 3; i++) {
    m_masses[i] = m_flavs[i].Mass();
  }
  // Neutrino: E_nu = (MW^2-ml^2)/(2*MW), cos(theta) in [-1,1], phi in [0,2*pi]
  double E_nu = (MW2-sqr(m_masses[1]))/(2.*MW);
  m_1to2moms[2] = E_nu*Vec4D(1.,0.,0.,-1.);
  // Electron: E_l = (MW^2+ml^2)/(2*MW), p_l = -p_nu
  double E_l = (MW2+sqr(m_masses[1]))/(2.*MW);
  m_1to2moms[1] = Vec4D(E_l,-Vec3D(m_1to2moms[2]));

  if ((m_1to2moms[1]+m_1to2moms[2])[0]-m_1to2moms[0][0] > 1E-6*m_1to2moms[0][0]) {
    msg_Out() << "You idiot. Momentum not conserved ... \nIncoming " << m_1to2moms[0] << "\nOutgoing " << m_1to2moms[1]+m_1to2moms[2];
  }
  

  // Max energy when u -> 0: 0 = Sqrt(M_W^2) - Sqrt(m_l^2) - Sqrt(Sqr(kappa_N)) - K0max.
  // kappa_N == K if neutrino takes all the recoil. Hence K0max = 0.5*(M_W - m_l)
  m_max = 0.5*(m_masses[0]-m_masses[1]);

  Complex res(0.,0.);

  // Boost particles into dipole CMS (p_W+p_l)
  Poincare DipCMS(m_1to2moms[0]+m_1to2moms[1]);
  Vec4D temp = m_1to2moms[0];
  DipCMS.Boost(temp);
  m_moms0[0] = temp;
  temp = m_1to2moms[1];
  DipCMS.Boost(temp);
  m_moms0[1] = temp;
  temp = m_1to2moms[2];
  DipCMS.Boost(temp);
  m_moms0[2] = temp;
  // Rotate such that the W lies in the -z-direction, lepton in z-direction  
  Poincare rot(m_moms0[0],Vec4D(0.,0.,0.,-1.));
  rot.Rotate(m_moms0[0]);
  rot.Rotate(m_moms0[1]);
  rot.Rotate(m_moms0[2]);


  bool passedcuts = true;
  double entry = 1.;
  double born = 0.;

  // Generate photon and correct other momenta
  // naming convention: m_moms0 are uncorrected momenta, m_moms are corrected momenta
  Setup_Corrected_Momenta(p);

  // Pseudorapidities and pTs have to be calculated in the lab frame
  // Boost into CMS of W after mapping, boost back into old dipole CMS frame,
  // rotate back and boost back into lab frame

  Poincare Ptop(m_moms[1]+m_moms[2]+m_moms[3]);
  Poincare Qtop(m_moms0[1]+m_moms0[2]);
  Vec4D Wmom = m_moms[0];
  Ptop.Boost(Wmom);
  Qtop.BoostBack(Wmom);
  rot.RotateBack(Wmom);
  DipCMS.BoostBack(Wmom);
  Vec4D leptonmom = m_moms[1];
  Ptop.Boost(leptonmom);
  Qtop.BoostBack(leptonmom);
  rot.RotateBack(leptonmom);
  DipCMS.BoostBack(leptonmom);  
  Vec4D neutrinomom = m_moms[2];
  Ptop.Boost(neutrinomom);
  Qtop.BoostBack(neutrinomom);
  rot.RotateBack(neutrinomom);
  DipCMS.BoostBack(neutrinomom);  
  Vec4D labphoton = m_moms[3];
  Ptop.Boost(labphoton);
  Qtop.BoostBack(labphoton);
  rot.RotateBack(labphoton);
  DipCMS.BoostBack(labphoton);

  Vec4D MomConservation = m_1to2moms[0]-(leptonmom+neutrinomom+labphoton);

  // Cut on lepton-photon invariant mass
  if ((leptonmom+labphoton).Abs() < 1.) {
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // lepton-neutrino invariant mass: < sqrt(s)
  // Invariant mass can be checked in any frame
  if ((leptonmom+neutrinomom).Abs() > m_1to2moms[0].Abs()) {
    msg_Out() << "Huh? Lepton invariant mass " << (leptonmom+neutrinomom).Abs() << " partonic COM energy " << m_1to2moms[0].Abs() << "\n";
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // lepton photon invariant mass: < sqrt(s)
  if (((leptonmom+neutrinomom+labphoton).Abs()-m_1to2moms[0].Abs()) > 1E-6*m_1to2moms[0].Abs()) {
    msg_Out() << "Huh? Lepton-photon invariant mass " << (leptonmom+neutrinomom+labphoton).Abs() << " partonic COM energy " << m_1to2moms[0].Abs() << "\n";
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // Check momentum conservation in lab frame
  if (MomConservation[0] > 1.E-6*m_1to2moms[0].Abs() ||
	   MomConservation[1] > 1.E-6*m_1to2moms[0].Abs() ||
	   MomConservation[2] > 1.E-6*m_1to2moms[0].Abs() ||
	   MomConservation[3] > 1.E-6*m_1to2moms[0].Abs()) {
    msg_Out() << "Momentum not conserved after boost back. Diff " << MomConservation << "\n";
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }
  // Check lepton momenta are larger than 0
  if (leptonmom[0] < 0. || neutrinomom[0] < 0.) {
    msg_Out() << "Lepton energy below 0.\nLepton " << leptonmom << "\nNeutrino "  << neutrinomom << std::endl;
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }
  // Check photon momentum is larger than 0
  if (labphoton[0] < 0.) {
    msg_Out() << "Photon energy below 0.\n " << labphoton << std::endl;;
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }

  if (passedcuts == false) {
    m_histos["E_lep"]->Insert(leptonmom[0],0.);
    m_histos["E_neu"]->Insert(neutrinomom[0],0.);
    m_histos["E_gamma"]->Insert(labphoton[0],0.);
    m_histos["Mlgamma"]->Insert((labphoton+leptonmom).Abs(),0.);
    m_histos["E_gamma_small"]->Insert(labphoton[0],0.);
    m_histos["Mlgamma_small"]->Insert((labphoton+leptonmom).Abs(),0.);
    m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(leptonmom))/(Vec3D(labphoton).Abs()*Vec3D(leptonmom).Abs()),0.);
  }
  

  
  if (passedcuts == true) {
      // Use generated and rescaled momenta in real
      m_1to3moms[0] = m_moms[0];
      m_1to3moms[1] = m_moms[1];
      m_1to3moms[2] = m_moms[2];
      m_1to3moms[3] = m_moms[3];


      
      Complex Decay[3];
      Complex LO_Decay[3];
      Complex sum_LO = Complex(0.,0.);
      Complex sum = Complex(0.,0.);
      double s = m_1to2moms[0].Abs2();
      
      for (unsigned int l=0; l<=number_pol; l++) {       // spin W in conj. amplitude
	for (unsigned int i=0; i<=1; i++) {	    // spin l.Bar
	  for (unsigned int j=0; j<=1; j++) {	    // spin l
	    m_1to3spins[2] = i;
	    m_1to3spins[1] = j;
	    for (unsigned int o=0; o<=1; o++) {   // spin G
	      Complex Dec = InfraredSubtractedME_1_05(l,o);
	      Complex conjDec = conj(InfraredSubtractedME_1_05(l,o));
	      Decay[l] += Dec*conjDec;
	    }
	  }
	}
	sum += 1./(16.*M_PI*M_PI*M_PI)*Decay[l];
      }

      for (unsigned int l=0; l<=number_pol; l++) {       // spin W in conj. amplitude
	for (unsigned int i=0; i<=1; i++) {	    // spin l.Bar
	  for (unsigned int j=0; j<=1; j++) {	    // spin l
	    m_1to2spins[2] = i;
	    m_1to2spins[1] = j;
	    Complex Dec = InfraredSubtractedME_0_0(l);
	    Complex conjDec = conj(InfraredSubtractedME_0_0(l));
	    LO_Decay[l] += Dec*conjDec;
	  }
	}
	sum_LO += LO_Decay[l];
      }


      double weight_jac_L = Weight_Jac_L();
      double weight_jac_M = Weight_Jac_M();
      double Uncorrected_Smod = Smod_Uncorr();

      // Safeguard against numerical issues ...
      if (Uncorrected_Smod <= pow(10.,-25.)  || m_angle_weight <= pow(10.,-25.)) {
	res = 0.;
      }
      else {
	res = sum/(3.*sum_LO)*(1.-sqr(m_masses[1])/(2.*MW2)-pow(m_masses[1],4.)/(2.*sqr(MW2)))*m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*weight_jac_L*weight_jac_M/Uncorrected_Smod/CouplingFactor(0,2);
	m_histos["E_lep"]->Insert(leptonmom[0],res.real());
	m_histos["E_neu"]->Insert(neutrinomom[0],res.real());
	m_histos["E_gamma"]->Insert(labphoton[0],res.real());
	m_histos["Mlgamma"]->Insert((labphoton+leptonmom).Abs(),res.real());
	m_histos["E_gamma_small"]->Insert(labphoton[0],res.real());
	m_histos["Mlgamma_small"]->Insert((labphoton+leptonmom).Abs(),res.real());
	m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(leptonmom))/(Vec3D(labphoton).Abs()*Vec3D(leptonmom).Abs()),res.real());
      }
      if (std::isinf(res.real()) || res.real() < 0. || std::isnan(res.real())) PRINT_VAR(res << " , " << log(m_max/m_omega) << " , " << m_angle_weight << " , " << GetBeta_1_1()/GetBeta_0_0() << " , " << weight_jac_L << " , " << weight_jac_M << " , " << Uncorrected_Smod);    
  }
  m_total_weight+=res;
  m_total_weight2+=sqr(res);
  return 1.;
}

void W_Decay_Real::Setup_Corrected_Momenta(const ATOOLS::Vec4D_Vector& p) {
  // Generate photon
  bool reject = true;
  double m = m_masses[1];
  double m2 = sqr(m);
  double M2 = sqr(m_masses[0]);
  while (reject) {
    // Generate photon
    m_moms[3] = Generate_One_Photon(m_moms0[0],m_moms0[1]);
    
    DetermineU(m_moms0[1],m_moms0[2]);
    // Ensure u is in [0,1], else regenerate
    if (m_u > 0. && m_u < 1.) reject = false;
  }
  // Rescale dipole momenta (taking into account masses correctly)
  Vec3D mom3 = m_u*Vec3D(m_moms0[0]);
  double momE = sqrt(M2+mom3.Sqr());
  m_moms[0] = Vec4D(momE,mom3);
  mom3 = m_u*Vec3D(m_moms0[1]);
  momE = sqrt(m2+mom3.Sqr());
  m_moms[1] = Vec4D(momE,mom3);
  // Compensate momenta with neutrino-momentum
  Vec3D kappaN = Vec3D(m_moms[3]);
  mom3 = m_u*Vec3D(m_moms0[2])-kappaN;
  momE = sqrt(mom3.Sqr());
  m_moms[2] = Vec4D(momE,mom3);

  // make sure W is still on-shell
  if (((m_moms[1]+m_moms[2]+m_moms[3]).Abs()-m_1to2moms[0].Abs()) > 1E-6*m_1to2moms[0].Abs()) {
    PRINT_VAR("\n" << m_u*Vec3D(m_moms0[0]) << " , " << sqrt(M2+(m_u*Vec3D(m_moms0[0])).Sqr()) << "\n" << m_u*Vec3D(m_moms0[1]) << " , " << sqrt(m2+(m_u*Vec3D(m_moms0[1])).Sqr()) << " , " << m_moms[1].Abs() << "\n" << m_u*Vec3D(m_moms0[2])-Vec3D(m_moms[3]) << " , " << sqrt((m_u*Vec3D(m_moms0[2])-Vec3D(m_moms[3])).Sqr()) << " , " << m_moms[2].Abs2() << "\n" << Vec3D(m_moms[3]) << " , " << m_moms[3].Abs2() << "\n" << m_moms[3][0] << "\n" << sqrt(Vec3D(m_moms[3]).Sqr()));
  }
  if (abs(m_moms[0].Abs()-m_moms0[0].Abs()) > 1.E-6*m_moms[0][0]) {
    msg_Out() << METHOD << "decay mass moved off-shell! before " << m_moms0[0].Abs() << " after " << m_moms[0].Abs() << " diff " << abs(m_moms0[0].Abs()-m_moms[0].Abs()) << " \n";
    msg_Out() << METHOD << "W Energy before " << m_moms0[0][0] << " after " << m_moms[0][0] << " photon momentum " << m_moms[3] << " u " << m_u << "\n";
  }
}

void W_Decay_Real::DetermineU(const Vec4D& QC, const Vec4D& QN) {
  // Use Newton-Raphson method
  bool closeenough = false;
  double m = m_masses[1];
  double m2 = sqr(m);
  double M2 = sqr(m_masses[0]);
  double f0,g0;
  double x = 0.5;
  double iter = 0;
  while (!closeenough) {
    // shortcuts for often used variables: QC2 = square of QC, QN2 = square of QN2, EW energy of W
    // El = energy of lepton, En = energy of neutrino
    double QC2 = Vec3D(QC)*Vec3D(QC);
    double QN2 = Vec3D(QN)*Vec3D(QN);
    double EW = sqrt(M2+sqr(x)*QC2);
    double El = sqrt(m2+sqr(x)*QC2);
    double En = sqrt(Vec3D(x*QN-m_moms[3])*Vec3D(x*QN-m_moms[3]));

    // For comparison with Mathematica - checked 24/09/2018
    msg_Debugging() << "Solve[Sqrt[x^2*" << QC2 << "`16+" << M2 << "`16]-Sqrt[x^2*" << QC2 << "`16+" << m2 << "`16]-Sqrt[x^2*" << QN2 << "`16-x*2*" << Vec3D(QN)*Vec3D(m_moms[3]) << "`16+" << Vec3D(m_moms[3])*Vec3D(m_moms[3]) << "`16]-" << m_moms[3][0] << "`16 == 0,x]" << std::endl;

    // Estimate of sum of E-components at x
    double f0 = EW - El - En - m_moms[3][0];
    // Consider result close to 0 if f0 is within 10^(-14) of 0
    if (dabs(f0) < pow(10.,-14.)*QC[0]) closeenough = true;
    if (closeenough == true && x < 0.) {
      closeenough = false;
      x = 0.8;
    }
    iter++;
    if (iter > 10) break;
    // Derivative of f0 at x
    g0 = x*QC2/EW - x*QC2/El - (x*QN2-Vec3D(m_moms[3])*Vec3D(QN))/En;
    // Update x
    if (g0 != 0.) x -= f0/g0;
    else {
      msg_Out() << METHOD << "Encountered f'(x) = 0.\n";
      x = -1.;
      break;
    }
  }
  if (x > 1. || x < 0.) {
    msg_Debugging() << "u = " << x << "    Returning m_u = -1." << std::endl;
    m_u = -1.;
    return;
  }
  if (x > 0. && x < 1.) {
    msg_Debugging() << "u = " << x << std::endl;
    m_u = x;
    return;
  }
}  

// Weight due to Lorentz trafo into Dipole CMS
double W_Decay_Real::Weight_Jac_L() {
  double QC0 = (m_moms0[1])[0];
  double PC0 = (m_moms[1])[0];
  // In dipole rest frame hence E = M
  double mMQ = (m_moms0[0]+m_moms0[1])[0];
  double mMP = (m_moms[0]+m_moms[1])[0];
  double QN0 = (m_moms0[2])[0];
  double PN0 = (m_moms[2])[0];
  double K0 = m_moms[3][0];
  return pow(mMP,3.)/pow(mMQ,3.)*(QC0+QN0)/(PC0+PN0+K0);
}

// Weight due to momentum mapping
double W_Decay_Real::Weight_Jac_M() {
  double sumq     = 0.;
  double sumpq     = 0.;
  double prod     = 1.;
  double N = 2;
  Vec3D Q         = Vec3D(0.,0.,0.);
  Vec3D P         = Vec3D(0.,0.,0.);
  for (int i(1); i < 3; i++) {
    sumq = sumq - Vec3D(m_moms0[i]).Sqr()/m_moms0[i][0];
    sumpq = sumpq - Vec3D(m_moms0[i])*Vec3D(m_moms[i])/m_moms[i][0];
    prod = prod*m_moms0[i][0]/m_moms[i][0];
    if (i == 1) {
      Q = Q + Vec3D(m_moms0[i]);
      P = P + Vec3D(m_moms[i]);
    }
  }
  sumq = sumq + (Q*Q)/sqrt(m_masses[0]*m_masses[0] + Q*Q);
  sumpq = sumpq + (P*Q)/sqrt(m_masses[0]*m_masses[0] + P*P);
  return pow(m_u,3.*N-4.) * sumq/sumpq * prod;
}


double W_Decay_Real::Generate_Energy(const double& omega, const double& max) {
  // Generate photon energy following 1/E
  return omega*pow(max/omega,ran->Get());
}

void W_Decay_Real::Generate_Photon_Angle(const Vec4D& p1, const Vec4D& p2) {
  Vec4D p1copy = p1;
  Vec4D p2copy = p2;
  // Need to switch order since p1 lies in -z-direction
  double b2 = Vec3D(p1copy).Abs()/p1copy[0];
  double b1 = Vec3D(p2copy).Abs()/p2copy[0];
  GenerateDipoleAngle(b1,b2);
  GenerateNullVector();
}

void W_Decay_Real::GenerateDipoleAngle(const double& b1, const double& b2) {
  // Generation of theta for two massive particles
  double P  = log((1.+b1)/(1.-b1))
                /(log((1.+b1)/(1.-b1))+log((1.+b2)/(1.-b2)));
  m_c = 0.;
  double a;
  if (ran->Get() < P) {
    double rnd = ran->Get();
    a   = 1./b1*log((1.+b1)/(1.-b1));
    m_c        = 1./b1*(1.-(1.+b1)*exp(-a*b1*rnd));
  }
  else {
    double rnd = ran->Get();
    a   = 1./b2*log((1.+b2)/(1.-b2));
    m_c        = 1./b2*((1.-b2)*exp(a*b2*rnd)-1.);
  }
  double weight = 1.-((1.-b1*b1)/((1.-b1*m_c)*(1.-b1*m_c))
  		      +(1.-b2*b2)/((1.+b2*m_c)*(1.+b2*m_c)))
    /(2.*(1.+b1*b2)/((1.-b1*m_c)*(1.+b2*m_c)));

  m_theta = acos(m_c);
  m_phi   = 2.*M_PI*ran->Get();

  // Combination of all weights due to photon angles - 2*M_PI due to
  // phi integration, part_weight for generation due to one term in
  // interference, weight as correction weight to full Eikonal, final
  // part is full integral of theta integration
  m_angle_weight = 2.*M_PI*weight*((1.+b1*b2)/(b1+b2)*log((1.+b1)*(1.+b2)/((1.-b1)*(1.-b2)))-2.);
  if (m_c < -1. || m_c > 1.) {
    msg_Out() << "Encountered c outside [-1,1], c = " << m_c << ". Theta = " << m_theta << "\n";
  }
}

void W_Decay_Real::GenerateNullVector() {
  // Setup null-vector in direction given by angles
  m_dir = Vec4D(1., sin(m_theta)*cos(m_phi), sin(m_theta)*sin(m_phi), cos(m_theta));
}

Vec4D W_Decay_Real::Generate_One_Photon(const Vec4D& p1, const Vec4D& p2) {
  double E = Generate_Energy(m_omega,m_max);
  if (E < m_below*m_deltas*p1.Abs()/2.) {
    msg_Out() << METHOD << ": This should not happen. Photon energy less than cut allows.\n";
    THROW(fatal_error, "cannot continue");
  }

  Generate_Photon_Angle(p1,p2);
  return E*m_dir;
}


Complex W_Decay_Real::InfraredSubtractedME_0_0() {
  Vec4C epsW = Polarization_Vector(m_moms0[0])[m_spins0[0]];
  XYZFunc XYZ(3,m_moms0,m_flavs,false);
  if (m_flavs[1].IsAnti()) {
    return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_spins0[1],epsW,2,m_spins0[2],0.,1.);
  }
  else {
    return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(2,m_spins0[2],epsW,1,m_spins0[1],0.,1.);
  }
}


Complex W_Decay_Real::InfraredSubtractedME_0_0(const int& Wpol) {
  Vec4C epsW = Polarization_Vector(m_moms0[0])[Wpol];
  XYZFunc XYZ(3,m_moms0,m_flavs,false);
  if (m_flavs[1].IsAnti()) {
    return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_1to2spins[1],epsW,2,m_1to2spins[2],0.,1.);
  }
  else {
    return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(2,m_1to2spins[2],epsW,1,m_1to2spins[1],0.,1.);
  }
}


Complex W_Decay_Real::GetBeta_0_0() {
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=number_pol; k++) {       // spin W reduce to 1 to check
        m_spins0[0] = k;
        m_spins0[1] = j;
        m_spins0[2] = i;
        Complex M_0_0 = InfraredSubtractedME_0_0();
        sum = sum + (M_0_0*conj(M_0_0))// .real()
	  ;
      }
    }
  }
  // spin avarage over initial state
  sum = (1./((double)number_pol+1.))*sum;
  return sum;
}

double W_Decay_Real::Smod() {
  // Smod calculated in case of one additional photon
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = m_flavs[0].Charge();
  double Zj = m_flavs[1].Charge();
  int    ti = -1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double W_Decay_Real::Smod_Vecs(const Vec4D& pW, const Vec4D& pl,const Vec4D& pnu,const Vec4D& pg) {
  // Smod calculated using input vectors - only used for tests
  double Zi = m_flavs[0].Charge();
  double Zj = m_flavs[1].Charge();
  int    ti = -1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pW/(pW*pg)-pl/(pl*pg)).Abs2();
}

double W_Decay_Real::Smod_Uncorr() {
  // Smod calculated in case of one additional photon using uncorrected momenta
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms0[0];
  Vec4D pj  = m_moms0[1];
  double Zi = m_flavs[0].Charge();
  double Zj = m_flavs[1].Charge();
  int    ti = -1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex W_Decay_Real::InfraredSubtractedME_1_05() {
  Vec4D real_moms[8];
  // Copy momenta for sanity ...
  for (size_t j = 0; j < 4; j++) real_moms[j] = m_moms[j];
  Vec4C epsW   = Polarization_Vector(real_moms[0])[m_spins[0]];
  Vec4C epsP   = conj(Polarization_Vector(real_moms[3])[m_spins[3]]);
  // Vec4C epsP   = m_moms[3]; // To check Ward identity
  Vec4D q      = real_moms[1]+real_moms[3];       // fermion propagator momenta
  double q2    = q.Abs2();
  Vec4D Q      = real_moms[0]-real_moms[3];       // boson propagator momenta
  double Q2    = Q.Abs2();
  double m     = m_masses[1];               // fermion mass/propagator pole
  double m2    = sqr(m);
  double M     = m_masses[0];               // boson mass/propagator pole
  double M2    = sqr(M);
  real_moms[4]    = real_moms[5] = q;             // enter those into m_moms
  m_flavs[4]   = m_flavs[1];                // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[1].Bar();
  XYZFunc XYZ(6,real_moms,m_flavs,false);
  // two diagrams
  // M_1 = -ie^2/(sqrt(2)sW) * 1/((pl+k)^2-m^2)
  //       * ubar(l)gamma^mu(-pl-k+m)gamma^nu P_L v(nu) eps_nu^W eps_mu^y*
  // M_2 = ie^2/(sqrt(2)sW) * 1/(pW-k)^2-M^2)
  //       * ubar(l)gamma_rho P_L v(nu)
  //       * [-2g^{rho,nu}pW^mu + g^{rho,mu}(pW-2k)^nu
  //          + 2g^{nu,mu}k^rho + 1/pW^2(pW-k)^rho pW^nu pW^mu]
  //       * eps_nu^W eps_mu^y*
  // Note that we need to use pW^2 instead of MW^2 in the propagator since W may be off-shell!
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Lorentz_Ten3C ten31,ten32,ten33,ten34,ten35;
  for (unsigned int s=0; s<=1; s++) {
    r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)
          *XYZ.X(4,s,epsW,2,m_spins[2],0.,1.);
    r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)
          *XYZ.X(5,s,epsW,2,m_spins[2],0.,1.);
  }
  Vec4D p = real_moms[0];
  Vec4D k = real_moms[3];
  // index ordering rho(1),nu(2),mu(3)
  // -2g^{rho,nu}pW^mu
  ten31 = BuildTensor(MetricTensor(),-2.*p);
  // g^{rho,mu}(pW-2.*k)^nu
  ten32 = BuildTensor(MetricTensor(),p-2.*k).Transpose(2,3);
  // 2g^{nu,mu}k^rho
  ten33 = BuildTensor(MetricTensor(),2.*k).Transpose(1,3);
  // 1/pW^2(pW-k)^rho pW^nu pW^mu
  ten34 = 1./(p*p)*BuildTensor(p-k,p,p);

  Lorentz_Ten3C ten = ten31+ten32+ten33+ten34;
  // v^\sigma = L^\sigma\mu\lambda epsW_\mu epsP_\lambda
  Vec4C v3 = Contraction(Contraction(ten,3,epsP),2,epsW);

  r3 = XYZ.X(1,m_spins[1],v3,2,m_spins[2],m_cR,m_cL);
  
  r1 *= (1.+m/sqrt(q2))/(2.*(q2-m2));
  r2 *= (1.-m/sqrt(q2))/(2.*(q2-m2));
  r3 *= -1./(Q2-p*p);

  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = Flavour(kf_none);
  // Note there was a factor of 2 wrong here - only applies to r1 and r2, since they involve spin sums
  return -(Complex(0.,1.)*m_e*m_e*CouplingFactor(0,1))/(sqrt(2.)*m_sW)*(r1+r2+r3);
}




Complex W_Decay_Real::InfraredSubtractedME_1_05(const int& Wpol, const int& Gpol) {
  Vec4D real_moms[6];
  // Copy momenta for sanity ...
  for (size_t j = 0; j < 4; j++) real_moms[j] = m_1to3moms[j];
  Vec4C epsW   = Polarization_Vector(real_moms[0])[Wpol];
  Vec4C epsP   = conj(Polarization_Vector(real_moms[3])[Gpol]);
  //Vec4C epsP   = real_moms[3]; // To check Ward identity
  Vec4D q      = real_moms[1]+real_moms[3];       // fermion propagator momenta
  double q2    = q.Abs2();
  Vec4D Q      = real_moms[0]-real_moms[3];       // boson propagator momenta
  double Q2    = Q.Abs2();
  double m     = m_masses[1];               // fermion mass/propagator pole
  double m2    = sqr(m);
  double M     = m_masses[0];               // boson mass/propagator pole
  double M2    = sqr(M);
  real_moms[4]    = real_moms[5] = q;             // enter those into m_moms
  m_real_flavs[4]   = m_real_flavs[1];                // set to corresponding particle/antiparticle
  m_real_flavs[5]   = m_real_flavs[1].Bar();
  XYZFunc XYZ(6,real_moms,m_real_flavs,false);
  // two diagrams
  // M_1 = -ie^2/(sqrt(2)sW) * 1/((pl+k)^2-m^2)
  //       * ubar(l)gamma^mu(-pl-k+m)gamma^nu P_L v(nu) eps_nu^W eps_mu^y*
  // M_2 = ie^2/(sqrt(2)sW) * 1/(pW-k)^2-M^2)
  //       * ubar(l)gamma_rho P_L v(nu)
  //       * [-2g^{rho,nu}pW^mu + g^{rho,mu}(pW-2k)^nu
  //          + 2g^{nu,mu}k^rho + 1/pW^2(pW-k)^rho pW^nu pW^mu]
  //       * eps_nu^W eps_mu^y*
  // Note that we need to use pW^2 instead of MW^2 in the propagator since W may be off-shell!
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Lorentz_Ten3C ten31,ten32,ten33,ten34,ten35;
  if (m_real_flavs[1].IsAnti()) {
    for (unsigned int s=0; s<=1; s++) {
      r1 += XYZ.X(1,m_1to3spins[1],epsP,4,s,1.,1.)
	*XYZ.X(4,s,epsW,2,m_1to3spins[2],0.,1.);
      r2 += XYZ.X(1,m_1to3spins[1],epsP,5,s,1.,1.)
	*XYZ.X(5,s,epsW,2,m_1to3spins[2],0.,1.);
    }
    Vec4D p = real_moms[0];
    Vec4D k = real_moms[3];
    // index ordering rho(1),nu(2),mu(3)
    // -2g^{rho,nu}pW^mu
    ten31 = BuildTensor(MetricTensor(),-2.*p);
    // g^{rho,mu}(pW-2.*k)^nu
    ten32 = BuildTensor(MetricTensor(),p-2.*k).Transpose(2,3);
    // 2g^{nu,mu}k^rho
    ten33 = BuildTensor(MetricTensor(),2.*k).Transpose(1,3);
    // 1/pW^2(pW-k)^rho pW^nu pW^mu
    ten34 = 1./(p*p)*BuildTensor(p-k,p,p);
    
    Lorentz_Ten3C ten = ten31+ten32+ten33+ten34;
    // v^\sigma = L^\sigma\mu\lambda epsW_\mu epsP_\lambda
    Vec4C v3 = Contraction(Contraction(ten,3,epsP),2,epsW);
    
    r3 = XYZ.X(1,m_1to3spins[1],v3,2,m_1to3spins[2],m_cR,m_cL);
    r3 *= -1./(Q2-p*p);
  }
  else {
    for (unsigned int s=0; s<=1; s++) {
      r1 += XYZ.X(4,s,epsP,1,m_1to3spins[1],1.,1.)
	*XYZ.X(2,m_1to3spins[2],epsW,4,s,0.,1.);
      r2 += XYZ.X(5,s,epsP,1,m_1to3spins[1],1.,1.)
	*XYZ.X(2,m_1to3spins[2],epsW,5,s,0.,1.);
    }
    Vec4D p = real_moms[0];
    Vec4D k = real_moms[5];
    // index ordering rho(1),nu(2),mu(3)
    // -2g^{rho,nu}pW^mu
    ten31 = BuildTensor(MetricTensor(),-2.*p);
    // g^{rho,mu}(pW-2.*k)^nu
    ten32 = BuildTensor(MetricTensor(),p-2.*k).Transpose(2,3);
    // 2g^{nu,mu}k^rho
    ten33 = BuildTensor(MetricTensor(),2.*k).Transpose(1,3);
    // 1/pW^2(pW-k)^rho pW^nu pW^mu
    ten34 = 1./(p*p)*BuildTensor(p-k,p,p);
    
    Lorentz_Ten3C ten = ten31+ten32+ten33+ten34;
    // v^\sigma = L^\sigma\mu\lambda epsW_\mu epsP_\lambda
    Vec4C v3 = Contraction(Contraction(ten,3,epsP),2,epsW);
    
    r3 = XYZ.X(2,m_1to3spins[2],v3,1,m_1to3spins[1],m_cR,m_cL);
    r3 *= -1./(Q2-p*p);
  }
  
  r1 *= (1.+m/sqrt(q2))/(2.*(q2-m2));
  r2 *= (1.-m/sqrt(q2))/(2.*(q2-m2));
  // erase intermediate entries from m_flavs
  m_real_flavs[4] = m_real_flavs[5] = Flavour(kf_none);
  // Note there was a factor of 2 wrong here - only applies to r1 and r2, since they involve spin sums
  return -(Complex(0.,1.)*m_e*m_e*CouplingFactor(0,1))/(sqrt(2.)*m_sW)*(r1+r2+r3);
}







Complex W_Decay_Real::GetBeta_1_1() {
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=number_pol; k++) {       // spin W reduced to check
        for (unsigned int l=0; l<=1; l++) {     // spin gamma
          m_spins[0] = k;
          m_spins[1] = j;
          m_spins[2] = i;
          m_spins[3] = l;
	  Complex M_1_05 = InfraredSubtractedME_1_05();
          sum = sum + (M_1_05*conj(M_1_05))// .real()
	    ;
        }
      }
    }
  }
  // spin avarage over initial state
  sum = (1./((double)number_pol+1.))*sum;
  sum = 1./(16.*M_PI*M_PI*M_PI)*sum;
  return sum;
}

DECLARE_TREEME2_GETTER(W_Decay_Real,
		       "W_Decay_Real")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,W_Decay_Real>::
operator()(const Process_Info &pi) const
{
  Default_Reader reader;
  if (reader.GetValue<int>("EXTRAXS_W_Decay_Real",0) != 1) return NULL;
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[2].IsLepton() && fl[3].IsLepton() && fl[3]!=fl[2].Bar() && fl[3].LeptonFamily()==fl[2].LeptonFamily() && 
      fl[0].IsQuark()  && fl[1]!=fl[0].Bar() && fl[0].IsUptype()!=fl[1].IsUptype() && fl[0].IsDowntype()!=fl[1].IsDowntype() &&
      (fl[0].IntCharge()+fl[1].IntCharge()==fl[2].IntCharge()+fl[3].IntCharge())) {
    if (pi.m_maxcpl[1]==3 && pi.m_mincpl[1]==3) {
      return new W_Decay_Real(pi,fl);
    }
  }
  return NULL;
}



















// double W_Decay_Real::operator()
// (const ATOOLS::Vec4D_Vector& p)
// {
//   // None of these ever happen
//   // if ((p[2]+p[3]).Abs() < 1.) msg_Error() << "Born cut not respected. m_ll " << (p[2]+p[3]).Abs() << std::endl;
//   // if (p[2].PPerp() < 25. ||  p[3].PPerp() < 25.) msg_Error() << "Born cut not respected. p[2].pT = " << p[2].PPerp() << " , p[3].pT = " << p[3].PPerp() << std::endl;
//   // if (abs(p[3].Eta()) > 2.5) msg_Error() << "Born cut not respected. p[2].eta = " << p[2].Eta() << " , p[3].eta = " << p[3].Eta() << std::endl;


//   // Energy cutoff in parton cms. Need to go below nominal value since generation in different frame
//   m_omega = m_below*m_deltas*(p[0]+p[1]).Abs()/2.;

//   m_masses[0] = (p[0]+p[1]).Abs();
//   for (int i(1); i < 3; i++) {
//     m_masses[i] = m_flavs[i].Mass();
//   }

//   // Max energy when u -> 0: 0 = Sqrt(M_W^2) - Sqrt(m_l^2) - Sqrt(Sqr(kappa_N)) - K0max.
//   // kappa_N == K if neutrino takes all the recoil. Hence K0max = 0.5*(M_W - m_l)
//   m_max = 0.5*(m_masses[0]-m_masses[1]);

//   Complex res(0.,0.);

//   // Boost particles into dipole CMS (p_W+p_l = p_q+p_qbar+p_l)
//   // Remember: p[2] is neutrino, p[3] lepton
//   Poincare DipCMS(p[0]+p[1]+p[3]);
//   Vec4D temp = p[0]+p[1];
//   DipCMS.Boost(temp);
//   m_moms0[0] = temp;
//   temp = p[3];
//   DipCMS.Boost(temp);
//   m_moms0[1] = temp;
//   temp = p[2];
//   DipCMS.Boost(temp);
//   m_moms0[2] = temp;
//   // Rotate such that the W lies in the -z-direction, lepton in z-direction
//   Poincare rot(m_moms0[0],Vec4D(0.,0.,0.,-1.));
//   rot.Rotate(m_moms0[0]);
//   rot.Rotate(m_moms0[1]);
//   rot.Rotate(m_moms0[2]);


//   bool passedcuts = true;
//   double entry = 1.;
//   double born = 0.;
//   Setup_Corrected_Momenta(p);


//   //if ((m_moms[0]-(m_moms[1]+m_moms[2]+m_moms[3]))[0] > pow(10.,-10.)) PRINT_VAR("\nBefore\n" << m_moms0[0] << "\n" << m_moms0[1] << "\n" << m_moms0[2]  << "\n" << m_moms0[1]+m_moms0[2]<< "\n\nAfter\n" << m_moms[0] << "\n" << m_moms[1] << "\n" << m_moms[2] << "\n" << m_moms[3] << "\n" << m_moms[1]+m_moms[2]+m_moms[3] << "\n" << m_moms[0]-(m_moms[1]+m_moms[2]+m_moms[3]));

//   // Vectors after boosts Ptop and Qtop
//   Vec4D rW, rl, rnu, rg;
//   // Vectors after rotation back
//   Vec4D tW, tl, tnu, tg;




//   // Pseudorapidities and pTs have to be calculated in the lab frame
//   // Boost into CMS of W after mapping, boost back into old dipole CMS frame,
//   // rotate back and boost back into lab frame

//   Poincare Ptop(m_moms[1]+m_moms[2]+m_moms[3]);
//   Poincare Qtop(m_moms0[1]+m_moms0[2]);
//   Vec4D Wmom = m_moms[0];
//   Ptop.Boost(Wmom);
//   Qtop.BoostBack(Wmom);
//   rW = Wmom;
//   rot.RotateBack(Wmom);
//   tW = Wmom;
//   DipCMS.BoostBack(Wmom);
//   Vec4D leptonmom = m_moms[1];
//   Ptop.Boost(leptonmom);
//   Qtop.BoostBack(leptonmom);
//   rl = leptonmom;
//   rot.RotateBack(leptonmom);
//   tl = leptonmom;
//   DipCMS.BoostBack(leptonmom);  
//   Vec4D neutrinomom = m_moms[2];
//   Ptop.Boost(neutrinomom);
//   Qtop.BoostBack(neutrinomom);
//   rnu = neutrinomom;
//   rot.RotateBack(neutrinomom);
//   tnu = neutrinomom;
//   DipCMS.BoostBack(neutrinomom);  
//   Vec4D labphoton = m_moms[3];
//   Ptop.Boost(labphoton);
//   Qtop.BoostBack(labphoton);
//   rg = labphoton;
//   rot.RotateBack(labphoton);
//   tg = labphoton;
//   DipCMS.BoostBack(labphoton);

//   // Need to apply photon cut in partonic CMS - boost labphoton there
//   Poincare PartonicCMS(p[0]+p[1]);
//   Vec4D labphoton_CMS = labphoton;
//   Vec4D leptonmom_CMS = leptonmom;
//   Vec4D neutrinomom_CMS = neutrinomom;
//   PartonicCMS.Boost(labphoton_CMS);
//   PartonicCMS.Boost(leptonmom_CMS);
//   PartonicCMS.Boost(neutrinomom_CMS);
//   Vec4D MomConservation = p[0]+p[1]-(leptonmom+neutrinomom+labphoton);

//   // Generated photon energy is NOT in partonic CMS but in frame boosted wrt partonic CMS by k
//   // hence need to boost to partonic CMS and check that it is above cut and below max energy photon
//   // can carry away!
//   if (labphoton_CMS[0] < m_deltas*(p[0]+p[1]).Abs()/2. || labphoton_CMS[0] > (Wmom.Abs2()-sqr(m_masses[1]))/(2.*Wmom.Abs())
//       ) {
//     entry = 0.;
//     res = 0.;
//     passedcuts = false;
//   }
//   // Cut on lepton-neutrino invariant mass
//   else if ((leptonmom+neutrinomom).Abs() < 1.) {
//     res = 0.;
//     entry = 0.;
//     passedcuts = false;
//   }
//   // lepton-neutrino invariant mass: < sqrt(s)
//   else if ((leptonmom+neutrinomom).Abs() > (p[0]+p[1]).Abs()) {
//     msg_Out() << "Huh? Lepton invariant mass " << (leptonmom+neutrinomom).Abs() << " partonic COM energy " << (p[0]+p[1]).Abs() << "\n";
//     // Invariant mass can be checked in any frame
//     res = 0.;
//     entry = 0.;
//     passedcuts = false;
//   }
//   // lepton photon invariant mass: < sqrt(s)
//   else if (((leptonmom+neutrinomom+labphoton).Abs()-(p[0]+p[1]).Abs()) > 1E-6*(p[0]+p[1]).Abs()) {
//     msg_Out() << "Huh? Lepton-photon invariant mass " << (leptonmom+neutrinomom+labphoton).Abs() << " partonic COM energy " << (p[0]+p[1]).Abs() << "\n";
//     res = 0.;
//     entry = 0.;
//     passedcuts = false;
//   }
//   // Check momentum conservation in lab frame
//   else if (MomConservation[0] > 1.E-6*(p[0]+p[1]).Abs() ||
// 	   MomConservation[1] > 1.E-6*(p[0]+p[1]).Abs() ||
// 	   MomConservation[2] > 1.E-6*(p[0]+p[1]).Abs() ||
// 	   MomConservation[3] > 1.E-6*(p[0]+p[1]).Abs()) {
//     msg_Out() << "Momentum not conserved after boost back. Diff " << MomConservation << "\n";
//     entry = 0.;
//     res = 0.;
//     passedcuts = false;
//   }
//   // Check lepton momenta are larger than 0
//   else if (leptonmom[0] < 0. || neutrinomom[0] < 0.) {
//     msg_Out() << "Lepton energy below 0.\nLepton " << leptonmom << "\nNeutrino "  << neutrinomom << std::endl;
//     entry = 0.;
//     res = 0.;
//     passedcuts = false;
//   }
//   // Check photon momentum is larger than 0
//   else if (labphoton[0] < 0.) {
//     msg_Out() << "Photon energy below 0.\n " << labphoton << std::endl;;
//     entry = 0.;
//     res = 0.;
//     passedcuts = false;
//   }
//   // Check other cuts: |eta(l)| < 2.5; pT(l) < 25 GeV; pT(nu) < 25 GeV
//   else if (abs(leptonmom.Eta()) > 2.5 ||
// 	   leptonmom.PPerp() < 25. ||
// 	   neutrinomom.PPerp() < 25.) {
//     res= 0.;
//     entry = 0.;
//     passedcuts = false;
//   }

//   // // Cuts on cos(theta) in CMS frame
//   // else if ((Vec3D(labphoton_CMS)*Vec3D(leptonmom))/(Vec3D(labphoton_CMS).Abs()*Vec3D(leptonmom).Abs()) < -0.95 || (Vec3D(labphoton_CMS)*Vec3D(leptonmom))/(Vec3D(labphoton_CMS).Abs()*Vec3D(leptonmom).Abs()) > 0.95) {
//   //   res = 0.;
//   //   entry = 0.;
//   //   passedcuts = false;
//   // }

//   if (passedcuts == false) {
//     m_histos["pT_lep"]->Insert(leptonmom.PPerp(),0.);
//     m_histos["pT_neu"]->Insert(neutrinomom.PPerp(),0.);
//     m_histos["Mlnu"]->Insert((leptonmom+neutrinomom).Abs(),0.);
//     m_histos["E_gamma"]->Insert(labphoton[0],0.);
//     m_histos["pT_gamma"]->Insert(labphoton.PPerp(),0.);
//     m_histos["Mlgamma"]->Insert((labphoton+leptonmom).Abs(),0.);
//     m_histos["E_gamma_small"]->Insert(labphoton[0],0.);
//     m_histos["pT_gamma_small"]->Insert(labphoton.PPerp(),0.);
//     m_histos["Mlgamma_small"]->Insert((labphoton+leptonmom).Abs(),0.);
//     m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(leptonmom))/(Vec3D(labphoton).Abs()*Vec3D(leptonmom).Abs()),0.);
//   }
  

  
//   if (passedcuts == true) {
//     double CKM;
//     // Determine CKM element
//     if ((m_quark_flavs[0].Kfcode() == kf_d && m_quark_flavs[1].Kfcode() == kf_u) ||
// 	(m_quark_flavs[0].Kfcode() == kf_u && m_quark_flavs[1].Kfcode() == kf_d)) CKM = 0.975;
//     else if ((m_quark_flavs[0].Kfcode() == kf_s && m_quark_flavs[1].Kfcode() == kf_u) || 
// 	     (m_quark_flavs[0].Kfcode() == kf_u && m_quark_flavs[1].Kfcode() == kf_s)) CKM = 0.222;
//     else if ((m_quark_flavs[0].Kfcode() == kf_d && m_quark_flavs[1].Kfcode() == kf_c) || 
// 	     (m_quark_flavs[0].Kfcode() == kf_c && m_quark_flavs[1].Kfcode() == kf_d)) CKM = 0.222;
//     else if ((m_quark_flavs[0].Kfcode() == kf_s && m_quark_flavs[1].Kfcode() == kf_c) || 
// 	     (m_quark_flavs[0].Kfcode() == kf_c && m_quark_flavs[1].Kfcode() == kf_s)) CKM = 0.975;
//     if (!pol_dep) {
//       double s = (p[0]+p[1]).Abs2();
//       double t(0.);
//       // Determine t. Need to check different cases due to ordering of quarks ...
//       if (m_quark_flavs[0].IsUptype() && m_flavs[2].Charge() != 0.) {
// 	// PRINT_VAR("02");
// 	// t=(p[0]-p[2]).Abs2();
// 	t = (p[1]-neutrinomom).Abs2();
//       }
//       else if (m_quark_flavs[0].IsUptype() && m_flavs[2].Charge() == 0.) {
// 	// PRINT_VAR("03");
// 	// t=(p[0]-p[3]).Abs2();
// 	t = (p[1]-neutrinomom).Abs2();
//       }
//       else if (m_quark_flavs[0].IsDowntype() && m_flavs[2].Charge() != 0.) {
// 	// PRINT_VAR("12");
// 	// t=(p[1]-p[2]).Abs2();
// 	t = (p[0]-neutrinomom).Abs2();
//       }
//       else if (m_quark_flavs[0].IsDowntype() && m_flavs[2].Charge() == 0.) {
// 	// t=(p[1]-p[3]).Abs2();
// 	t = (p[0]-neutrinomom).Abs2();
//       }
      
//       // Born amplitude
//       chi = sqr(t)/(sqr(s-MW2) + GW2*MW2);
//       born = sqr(CKM)/(4.*abs(m_sW2)*abs(m_sW2))*chi*(1.-sqr(m_masses[1])/(2.*MW2)-pow(m_masses[1],4.)/(2.*sqr(MW2)));
      
//       // Calculate weights - jac_L: weight due to LT into dipole CMS,
//       // jac_M: weight due to momentum mapping, Uncorrected_Smod:
//       // eikonal factor using uncorrected lepton momenta
//       double weight_jac_L = Weight_Jac_L();
//       double weight_jac_M = Weight_Jac_M();
//       double Uncorrected_Smod = Smod_Uncorr();
      
//       // Safeguard against numerical issues ...
//       if (Uncorrected_Smod <= pow(10.,-25.)  || m_angle_weight <= pow(10.,-25.)) {
// 	res = 0.;
//       }
//       else {
// 	// Combine all weights: log(m_max/m_omega) due to energy
// 	// generation, m_angle_weight due to angle generation,
// 	// GetBeta_1_1()/GetBeta_0_0() due to ME, weight_jac_L,
// 	// weight_jac_M Jacobians, 1/Uncorrected_Smod to correct for
// 	// generation following "wrong" momenta
	
// 	// Vec4D beam_1 = Vec4D(   425.24169922997396      ,   0.0000000000000000      ,   0.0000000000000000      ,   425.24169922997396      );
// 	// Vec4D beam_2 = Vec4D(   3.5988798021605533      ,   0.0000000000000000      ,   0.0000000000000000      ,  -3.5988798021605533      );
// 	// m_moms[1] = Vec4D(   65.293073295119186      ,  -4.3239670770541236      ,  -25.644513345120743      ,   59.890296910844910      );
// 	// m_moms[2] = Vec4D(   356.11589515703719      ,   4.8161621112877215      ,   28.563377950023245      ,   354.93586576464725      );
// 	// m_moms[3] = Vec4D(   7.4316105799782033      , -0.49219503423359773      ,  -2.9188646049025011      ,   6.8166567523212649      );
// 	// m_moms[0] = m_moms[1]+m_moms[2]+m_moms[3];
// 	// m_moms0[0] = beam_1+beam_2;
// 	// m_moms0[1] = beam_1;
// 	// m_moms0[2] = beam_2;
// 	// s = (beam_1+beam_2).Abs2();
// 	// t = (beam_1-m_moms[2]).Abs2();
// 	// m_masses[0] = s;

// 	// PRINT_VAR(16.*pow(M_PI,3.)*sqr(t)/(sqr(s-MW2) + GW2*MW2)*GetBeta_1_1()/GetBeta_0_0()/(4.*abs(m_sW2)*abs(m_sW2))*sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.);
// 	// PRINT_VAR(16.*pow(M_PI,3.)*sqr(t)/(sqr(s-MW2) + GW2*MW2)*GetBeta_1_1()/GetBeta_0_0()/(4.*abs(m_sW2)*abs(m_sW2))*sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3./(64.*sqr(M_PI*m_alpha*CouplingFactor(0,1)/m_sW2)*sqr(M_PI)/3.));
// 	// PRINT_VAR((64.*sqr(M_PI*m_alpha*CouplingFactor(0,1)/m_sW2)*sqr(M_PI)/3.));

// 	// Think about making Beta_1_1/Beta_0_0 polarization dependent
// 	// res = sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*GetBeta_1_1()/GetBeta_0_0()*born*weight_jac_L*weight_jac_M/Uncorrected_Smod;

// 	// Phase space only
// 	res = 3./(16.*M_PI*M_PI*M_PI)*sqr(CKM)*m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*weight_jac_L*weight_jac_M/Uncorrected_Smod;
//       }
//     }
//     else {
//       m_quarkmoms[0] = p[0];
//       m_quarkmoms[1] = p[1];
//       DipCMS.Boost(m_quarkmoms[0]);
//       DipCMS.Boost(m_quarkmoms[1]);
//       rot.Rotate(m_quarkmoms[0]);
//       rot.Rotate(m_quarkmoms[1]);
//       m_2to2moms[0] = m_moms0[0];
//       if (IsEqual(m_quarkmoms[0].Abs(),m_born_flavs[1].Mass()),1e-6) {
// 	m_2to2moms[1] = m_quarkmoms[0];
// 	m_2to2moms[2] = m_quarkmoms[1];
//       }
//       else {
// 	m_2to2moms[1] = m_quarkmoms[1];
// 	m_2to2moms[2] = m_quarkmoms[0];
//       }
//       if (IsEqual(m_moms0[1].Abs(),m_born_flavs[4].Mass()),1e-6) { 
// 	m_2to2moms[3] = m_moms0[2];
// 	m_2to2moms[4] = m_moms0[1];
//       }
//       else {
// 	m_2to2moms[3] = m_moms0[1];
// 	m_2to2moms[4] = m_moms0[2];
//       }   

//       // // Use lab momenta in real
//       // m_2to3moms[0] = rW;
//       // if (IsEqual(p[0].Abs(),m_real_flavs[1].Mass()),1e-6) {
//       // 	m_2to3moms[1] = m_quarkmoms[0];
//       // 	m_2to3moms[2] = m_quarkmoms[1];
//       // }
//       // else {
//       // 	m_2to3moms[1] = m_quarkmoms[1];
//       // 	m_2to3moms[2] = m_quarkmoms[0];
//       // }
//       // if (IsEqual(leptonmom.Abs(),m_born_flavs[4].Mass()),1e-5) { 
//       // 	m_2to3moms[3] = rnu;
//       // 	m_2to3moms[4] = rl;
//       // }
//       // else {
//       // 	m_2to3moms[3] = rl;
//       // 	m_2to3moms[4] = rnu;
//       // }   
//       // m_2to3moms[5] = rg;
//       // Use lab momenta in real
//       m_2to3moms[0] = m_moms[0];
//       if (IsEqual(p[0].Abs(),m_real_flavs[1].Mass()),1e-6) {
// 	m_2to3moms[1] = m_quarkmoms[0];
// 	m_2to3moms[2] = m_quarkmoms[1];
//       }
//       else {
// 	m_2to3moms[1] = m_quarkmoms[1];
// 	m_2to3moms[2] = m_quarkmoms[0];
//       }
//       Qtop.Boost(m_2to3moms[1]);
//       Ptop.BoostBack(m_2to3moms[1]);
//       Qtop.Boost(m_2to3moms[2]);
//       Ptop.BoostBack(m_2to3moms[2]);
//       if (IsEqual(leptonmom.Abs(),m_born_flavs[4].Mass()),1e-5) { 
// 	m_2to3moms[3] = m_moms[2];
// 	m_2to3moms[4] = m_moms[1];
//       }
//       else {
// 	m_2to3moms[3] = m_moms[1];
// 	m_2to3moms[4] = m_moms[2];
//       }   
//       m_2to3moms[5] = m_moms[3];


      

//       // Spin sum test: Sum over W polarizations of ProdME()*DecayME() == Z(1,2,3,4)
//       // for (unsigned int r=0; r<=1; r++) {           // spin quark1
//       // 	for (unsigned int q=0; q<=1; q++) {	    // spin quark2
//       // 	  for (unsigned int i=0; i<=1; i++) {	    // spin l.Bar
//       // 	    for (unsigned int j=0; j<=1; j++) {	    // spin l
//       // 	      Complex sum = Complex(0.,0.);
//       // 	      for (unsigned int k=0; k<=number_pol; k++) {       // spin W reduce to 1 to check
//       // 		m_2to2spins[0] = k;
//       // 		m_2to2spins[1] = r;
//       // 		m_2to2spins[2] = q;
//       // 		m_2to2spins[3] = i;
//       // 		m_2to2spins[4] = j;
//       // 		Complex Pro = ProdME();
//       // 		Complex Dec = LO_DecayME();
//       // 		sum += Pro*Dec;
//       // 	      }
//       // 	      XYZFunc XYZ(5,m_2to2moms,m_born_flavs,false,1);
//       // 	      if (m_born_flavs[1].IsAnti() && m_born_flavs[3].IsAnti()) {
//       // 		PRINT_VAR(r << " , " << q << " , " << i << " , " << j);
//       // 		Complex Z = XYZ.Z(1,r,2,q,3,i,4,j,0.,1.,0.,1.)*sqr(m_e)*CouplingFactor(0,1)/(2.*m_sW2);
//       // 		PRINT_VAR("\nsum " << sum << "\nZ   " << Z);
//       // 	      }
//       // 	      else if (!m_born_flavs[1].IsAnti() && m_born_flavs[3].IsAnti()) {
//       // 		PRINT_VAR(q << " , " << r << " , " << i << " , " << j);
//       // 		Complex Z = XYZ.Z(2,q,1,r,3,i,4,j,0.,1.,0.,1.)*sqr(m_e)*CouplingFactor(0,1)/(2.*m_sW2);
//       // 		PRINT_VAR("\nsum " << sum << "\nZ   " << Z);
//       // 	      }
//       // 	      else if (m_born_flavs[1].IsAnti() && !m_born_flavs[3].IsAnti()) {
//       // 		PRINT_VAR(r << " , " << q << " , " << j << " , " << i);
//       // 		Complex Z = XYZ.Z(1,r,2,q,4,j,3,i,0.,1.,0.,1.)*sqr(m_e)*CouplingFactor(0,1)/(2.*m_sW2);
//       // 		PRINT_VAR("\nsum " << sum << "\nZ   " << Z);
//       // 	      }
//       // 	      else if (!m_born_flavs[1].IsAnti() && !m_born_flavs[3].IsAnti()) {
//       // 		PRINT_VAR(q << " , " << r << " , " << j << " , " << i);
//       // 		Complex Z = XYZ.Z(2,q,1,r,4,j,3,i,0.,1.,0.,1.)*sqr(m_e)*CouplingFactor(0,1)/(2.*m_sW2);
//       // 		PRINT_VAR("\nsum " << sum << "\nZ   " << Z);
//       // 	      }
//       // 	    }
//       // 	  }
//       // 	}
//       // }
    
//  // m_2to2moms[1] = Vec4D(   139.35120657629000      ,   0.0000000000000000      ,   0.0000000000000000      ,   139.35120657629000      );
//  // m_2to2moms[2] = Vec4D(   23.702209738518064      ,   0.0000000000000000      ,   0.0000000000000000      ,  -23.702209738518064      );
//  // m_2to2moms[4] = Vec4D(   58.505930776091034      ,  -17.357185765462447      ,   49.781128005871508      ,   25.367525154849034      );
//  // m_2to2moms[3] = Vec4D(   104.54748553871703      ,   17.357185765462447      ,  -49.781128005871508      ,   90.281471682922898      );



//       // PRINT_VAR("\n  " << (m_born_flavs[1].IsAnti()?"-":"") << m_born_flavs[1].Kfcode() << "  " << m_2to2moms[1][0] << "  " << m_2to2moms[1][1] <<  "  " << m_2to2moms[1][2] << "  " << m_2to2moms[1][3] << "\n  " << (m_born_flavs[2].IsAnti()?"-":"") << m_born_flavs[2].Kfcode() << "  " << m_2to2moms[2][0] << "  " << m_2to2moms[2][1] <<  "  " << m_2to2moms[2][2] << "  " << m_2to2moms[2][3] << "\n  " << (m_born_flavs[3].IsAnti()?"-":"") << m_born_flavs[3].Kfcode() << "  " << m_2to2moms[3][0] << "  " << m_2to2moms[3][1] <<  "  " << m_2to2moms[3][2] << "  " << m_2to2moms[3][3] << "\n  " << (m_born_flavs[4].IsAnti()?"-":"") << m_born_flavs[4].Kfcode() << "  " << m_2to2moms[4][0] << "  " << m_2to2moms[4][1] <<  "  " << m_2to2moms[4][2] << "  " << m_2to2moms[4][3]);


//       // m_2to2moms[0] = m_2to2moms[3]+m_2to2moms[4];
//       // Poincare Twoto2DipCMS = Poincare(m_2to2moms[0]+m_2to2moms[4]);
//       // Twoto2DipCMS.Boost(m_2to2moms[0]);
//       // Twoto2DipCMS.Boost(m_2to2moms[1]);
//       // Twoto2DipCMS.Boost(m_2to2moms[2]);
//       // Twoto2DipCMS.Boost(m_2to2moms[3]);
//       // Twoto2DipCMS.Boost(m_2to2moms[4]);
//       // m_masses[0] = m_2to2moms[0].Abs();
//       Complex Production[3][3];
//       Complex Decay[3][3];
//       Complex sum = Complex(0.,0.);
//       double s = (p[0]+p[1]).Abs2();
//       // double s = (m_2to2moms[1]+m_2to2moms[2]).Abs2();
      
//       // // LO setup
//       // for (unsigned int k=0; k<=number_pol; k++) {       // spin W 
//       // 	for (unsigned int l=0; l<=number_pol; l++) {       // spin W in conj. amplitude
//       // 	  for (unsigned int r=0; r<=1; r++) {           // spin quark1
//       // 	    for (unsigned int q=0; q<=1; q++) {	    // spin quark2
//       // 	      m_2to2spins[1] = r;
//       // 	      m_2to2spins[2] = q;
//       // 	      Complex Pro = ProdME(k);
//       // 	      Complex conjPro = conj(ProdME(l));//ProdMEconj(l);
//       // 	      Production[k][l] += Pro*conjPro;
//       // 	    }
//       // 	  }
//       // 	  for (unsigned int i=0; i<=1; i++) {	    // spin l.Bar
//       // 	    for (unsigned int j=0; j<=1; j++) {	    // spin l
//       // 	      m_2to2spins[3] = i;
//       // 	      m_2to2spins[4] = j;
//       // 	      Complex Dec = LO_DecayME(k);
//       // 	      Complex conjDec = conj(LO_DecayME(l));//LO_DecayMEconj(l);
//       // 	      Decay[k][l] += Dec*conjDec;
//       // 	    }
//       // 	  }
//       // 	  sum += Production[k][l]*Decay[k][l];
//       // 	}
//       // }

//       // res = sqr(CKM)*1./(sqr(s-MW2)+MW2*GW2)*sum/4.*(1.-sqr(m_masses[1])/(2.*sqr(m_masses[0]))-pow(m_masses[1],4.)/(2.*sqr(sqr(m_masses[0]))));
//       // // PRINT_VAR("\nres " << res/3. <<
//       // // 		"\nden " << sqr(s-MW2)+MW2*GW2 << 
//       // // 		"\no_fac " << sqr(M_PI*m_alpha*CouplingFactor(0,1)/m_sW2)/(sqr(s-MW2)+MW2*GW2)/3. <<
//       // // 		"\npsq " << m_2to2moms[4].Abs2());
//       // m_histos["pT_lep"]->Insert(p[3].PPerp(),res.real());
//       // m_histos["pT_neu"]->Insert(p[2].PPerp(),res.real());
//       // m_histos["Mlnu"]->Insert((p[2]+p[3]).Abs(),res.real());
      
//       for (unsigned int k=0; k<=number_pol; k++) {       // spin W 
//       	for (unsigned int l=0; l<=number_pol; l++) {       // spin W in conj. amplitude
//       	  for (unsigned int r=0; r<=1; r++) {           // spin quark1
//       	    for (unsigned int q=0; q<=1; q++) {	    // spin quark2
//       	      m_2to3spins[1] = r;
//       	      m_2to3spins[2] = q;
//       	      Complex Pro = RealProdME(k);
//       	      Complex conjPro = conj(RealProdME(l));//ProdMEconj(l);
//       	      Production[k][l] += Pro*conjPro;
//       	    }
//       	  }
//       	  for (unsigned int i=0; i<=1; i++) {	    // spin l.Bar
//       	    for (unsigned int j=0; j<=1; j++) {	    // spin l
//       	      m_2to3spins[3] = i;
//       	      m_2to3spins[4] = j;
//       	      for (unsigned int o=0; o<=1; o++) {   // spin G
//       		Complex Dec = InfraredSubtractedME_1_05(k,o);
//       		Complex conjDec = conj(InfraredSubtractedME_1_05(l,o));//LO_DecayMEconj(l);
//       		Decay[k][l] += Dec*conjDec;
//       	      }
//       	    }
//       	  }
//       	  sum += 1./(16.*M_PI*M_PI*M_PI)*Production[k][l]*Decay[k][l];
//       	}
//       }





//       // msg_Out() << "Production\n(" << Production[0][0] << " , " << Production[0][1] << " , " << Production[0][2] << ")\n(" << "(" << Production[1][0] << " , " << Production[1][1] << " , " << Production[1][2] << ")\n(" << "(" << Production[2][0] << " , " << Production[2][1] << " , " << Production[2][2] << ")\n\nDecay\n(" << Decay[0][0] << " , " << Decay[0][1] << " , " << Decay[0][2] << ")\n(" << "(" << Decay[1][0] << " , " << Decay[1][1] << " , " << Decay[1][2] << ")\n(" << "(" << Decay[2][0] << " , " << Decay[2][1] << " , " << Decay[2][2] << ")";

//       double weight_jac_L = Weight_Jac_L();
//       double weight_jac_M = Weight_Jac_M();
//       double Uncorrected_Smod = Smod_Uncorr();

//       // Safeguard against numerical issues ...
//       if (Uncorrected_Smod <= pow(10.,-25.)  || m_angle_weight <= pow(10.,-25.)) {
// 	// res = 0.;
//       }
//       else {
// 	res = sqr(CKM)*1./(sqr(s-MW2)+MW2*GW2)*sum/4.*(1.-sqr(m_masses[1])/(2.*MW2)-pow(m_masses[1],4.)/(2.*sqr(MW2)))*m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*weight_jac_L*weight_jac_M/Uncorrected_Smod;
// 	m_histos["pT_lep"]->Insert(leptonmom.PPerp(),res.real());
// 	m_histos["pT_neu"]->Insert(neutrinomom.PPerp(),res.real());
// 	m_histos["Mlnu"]->Insert((leptonmom+neutrinomom).MPerp(),res.real());
// 	m_histos["E_gamma"]->Insert(labphoton[0],res.real());
// 	m_histos["pT_gamma"]->Insert(labphoton.PPerp(),res.real());
// 	m_histos["Mlgamma"]->Insert((labphoton+leptonmom).Abs(),res.real());
// 	m_histos["E_gamma_small"]->Insert(labphoton[0],res.real());
// 	m_histos["pT_gamma_small"]->Insert(labphoton.PPerp(),res.real());
// 	m_histos["Mlgamma_small"]->Insert((labphoton+leptonmom).Abs(),res.real());
// 	m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(leptonmom))/(Vec3D(labphoton).Abs()*Vec3D(leptonmom).Abs()),res.real());
// 	// Phase space only 
// 	// res = 1./(16.*M_PI*M_PI*M_PI)*sqr(CKM)*m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*weight_jac_L*weight_jac_M/Uncorrected_Smod;
//       }
//     }
//     // m_histos[string("m_u")]->Insert(m_u,entry);
//     // m_histos["E_gamma_acc"]->Insert(m_moms[3][0],entry);
//     // m_histos["E_labgamma_acc"]->Insert(labphoton[0],entry);
//     // if (std::isinf(res.real()) || res.real() < 0. || std::isnan(res.real())) PRINT_VAR(res << " , " << log(m_max/m_omega) << " , " << m_angle_weight << " , " << GetBeta_1_1()/GetBeta_0_0() << " , " << weight_jac_L << " , " << weight_jac_M << " , " << Uncorrected_Smod);    
//   }

//   // 1/3 factor due to quark colour averaging
//   // return sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*born;
//   return /*sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*/1./3.*res.real();
// }
// Not currently used
// double W_Decay_Real::DetermineMaxEnergy() {
//   double sum(0.);
//   double m = m_masses[1];
//   double m2 = sqr(m);
//   double M2 = sqr(m_masses[0]);
//   sum += m;
//   unsigned int count = 0;
//   double fac         = 1./3.;
//   double x0          = 0.5*(m_masses[0]-sum);
//   double xNminus1    = x0;
//   double xN          = 0.;
//   double F_xNminus1  = 0.;
//   while (abs(xN - xNminus1) > 1E-6) {
//     if (count==500) {
//       msg_Out()<<"failed to determine maximum photon energy... set to IR cut-off..."<<endl;
//       return m_omega;
//     }
//     xNminus1 = xN;
//     F_xNminus1 = sqr(fac*xNminus1);
//     sum = 0.;
//     sum += sqrt(m2+F_xNminus1) + sqrt(F_xNminus1);
//     xN = sqrt(m_masses[0]*m_masses[0] + F_xNminus1) - sum;
//     count++;
//   }
//   if (xN<0.) return m_omega;
//   return xN;
// }






// Complex W_Decay_Real::ProdME(const int& Wpol) {
//   Vec4C epsW = conj(Polarization_Vector(m_2to2moms[0])[Wpol]);
//   XYZFunc XYZ(5,m_2to2moms,m_born_flavs,false,1);
//   if (m_born_flavs[1].IsAnti()) {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_2to2spins[1],epsW,2,m_2to2spins[2],0.,1.);
//   }
//   else {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(2,m_2to2spins[2],epsW,1,m_2to2spins[1],0.,1.);
//   }
// }

// Complex W_Decay_Real::RealProdME(const int& Wpol) {
//   Vec4C epsW = conj(Polarization_Vector(m_2to3moms[0])[Wpol]);
//   XYZFunc XYZ(5,m_2to3moms,m_real_flavs,false,1);
//   if (m_born_flavs[1].IsAnti()) {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_2to3spins[1],epsW,2,m_2to3spins[2],0.,1.);
//   }
//   else {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(2,m_2to3spins[2],epsW,1,m_2to3spins[1],0.,1.);
//   }
// }


// Complex W_Decay_Real::LO_DecayME(const int& Wpol) {
//   Vec4C epsW = Polarization_Vector(m_2to2moms[0])[Wpol];
//   XYZFunc XYZ(5,m_2to2moms,m_born_flavs,false,1);
//   if (m_born_flavs[4].IsAnti()) {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(3,m_2to2spins[3],epsW,4,m_2to2spins[4],0.,1.);
//   }
//   else {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(4,m_2to2spins[4],epsW,3,m_2to2spins[3],0.,1.);
//   }
// }



// Complex W_Decay_Real::ProdMEconj(const int& Wpol) {
//   Vec4C epsW = (Polarization_Vector(m_2to2moms[0])[Wpol]);
//   XYZFunc XYZ(5,m_2to2moms,m_born_flavs,false,1);
//   if (m_born_flavs[1].IsAnti()) {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(2,m_2to2spins[2],epsW,1,m_2to2spins[1],0.,1.);
//   }
//   else {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(1,m_2to2spins[1],epsW,2,m_2to2spins[2],0.,1.);
//   }
// }


// Complex W_Decay_Real::LO_DecayMEconj(const int& Wpol) {
//   Vec4C epsW = conj(Polarization_Vector(m_2to2moms[0])[Wpol]);
//   XYZFunc XYZ(5,m_2to2moms,m_born_flavs,false,1);
//   if (m_born_flavs[3].IsAnti()) {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(4,m_2to2spins[4],epsW,3,m_2to2spins[3],0.,1.);
//   }
//   else {
//     return Complex(0.,1.)*m_e*sqrt(CouplingFactor(0,1))/(sqrt(2.)*m_sW)*XYZ.X(3,m_2to2spins[3],epsW,4,m_2to2spins[4],0.,1.);
//   }
// }



