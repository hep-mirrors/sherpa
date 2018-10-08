#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
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

  class Z_Decay_Real : public ME2_Base {
  private:
    int                m_spins[4], m_test, m_counter;
    int                m_spins0[3];
    int                MC_method, m_soft_approx;
    Flavour    m_flavs[9], m_quark, m_lepton;
    double m_e, m_alpha, m_deltas, m_omega, m_mu2, m_below, m_prev;
    double chi1,chi2,qedterm,intterm,Zterm,m_kappa,qf,qe,vf,af,ve,ae,sin2tw,mass;
    Complex             m_sW, m_sW2, m_sW4;
    Complex             m_cW, m_cW2, m_cW4;
    double m_theta, m_phi, m_c, m_max, m_u;
    double MW2;
    double MZ2, GZ2;
    double m_angle_weight;
    Complex muW2;
    Complex muZ2;
    ATOOLS::Vec4D m_dir;
    ATOOLS::Vec4D m_moms[4];
    ATOOLS::Vec4D m_moms0[3];
    double m_masses[3];
    Complex    m_cL[2];
    Complex    m_cR[2];
    // EXTRAXS::ME2_Base * p_bornme;
    std::map<std::string, ATOOLS::Histogram * >m_histos;
    std::map<std::string, ATOOLS::Histogram_2D * >m_histos_2D;
  public:

    Z_Decay_Real(const Process_Info& pi,const Flavour_Vector& fl);
    ~Z_Decay_Real();
    Complex InfraredSubtractedME_0_0(const int& b);
    Complex InfraredSubtractedME_1_05(const int& b);
    Complex GetBeta_0_0(const int& b1, const int& b2);
    Complex GetBeta_1_1(const int& b1, const int& b2);
    double Smod();
    double Smod_Uncorr();
    double Smod_Norm();
    ATOOLS::Vec4D Generate_One_Photon(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    double Generate_Energy(const double& omega, const double& max);
    void Generate_Photon_Angle(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    void GenerateDipoleAngle(const double& b1, const double& b2);
    void GenerateNullVector();
    double Weight_Jac_M();
    double Weight_Jac_L();
    double operator()(const ATOOLS::Vec4D_Vector& mom);
    void InitAnalysis();
    void FinishAnalysis(const Flavour& fl, const Flavour& fl_lep);
    void Setup_Corrected_Momenta(const ATOOLS::Vec4D_Vector& mom);
    double PDF(const double& c, const double& b);
    double CDF(const double& c, const double& b);
  };
}

Z_Decay_Real::Z_Decay_Real
(const Process_Info& pi, const Flavour_Vector& fl) : ME2_Base(pi, fl)
{
  Data_Reader reader(" ",";","#","=");
  // Read in cutoff
  m_deltas = reader.GetValue<double>("YFS_DELTAS",0.001);
  string widthscheme=reader.GetValue<string>("WIDTH_SCHEME","CMS");
  MC_method = reader.GetValue<int>("EXTRAXS_MC_METHOD",1);
  m_soft_approx = reader.GetValue<int>("EXTRAXS_SOFT_ONLY",0);
  Complex I = Complex(0.,1.);
  // Set up masses, couplings for born calculation
  sin2tw = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  mass     = fl[2].Mass();
  qe	   = fl[0].Charge();
  qf	   = fl[2].Charge();
  ae	   = fl[0].IsoWeak();	   
  af	   = fl[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  m_kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  
  m_quark = fl[0];
  m_lepton = fl[2];

  double  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  MW2 = pow(MW,2.);
  MZ2 = pow(MZ,2.);
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  GZ2 = GZ*GZ;
  muW2 = MW*(MW-I*GW);
  muZ2 = MZ*(MZ-I*GZ);
  if (widthscheme == "CMS") {
    m_cW2=muW2/muZ2;
    m_sW2=1.-m_cW2;
  }
  else if (widthscheme == "Fixed") {
    // This is the scheme in which all calculations are performed
    m_cW2=MW2/MZ2;
    m_sW2=1.-m_cW2;
  }
  m_sW = sqrt(m_sW2);
  m_cW = sqrt(m_cW2);
  m_e = sqrt(4.*M_PI*AlphaQED());
  m_alpha = AlphaQED();
  // Set up flavours for decay
  m_flavs[0]  = Flavour(kf_Z);
  if (!fl[2].IsAnti()) {
    m_flavs[1] = fl[2]; 
    m_flavs[2] = fl[3]; 
  }
  else {
    m_flavs[2] = fl[2]; 
    m_flavs[1] = fl[3]; 
  }
  m_flavs[3] = Flavour(kf_photon);
  // left/right-handed couplings to be used in X-functions - [0] for
  // couplings to photon, [1] for couplings to Z
  m_cL[0] = -I*m_e*m_flavs[1].Charge();
  m_cR[0] = -I*m_e*m_flavs[1].Charge();
  m_cL[1] = I*m_e/(2.*m_sW*m_cW)*(2.*m_flavs[1].IsoWeak()
				 -2.*m_flavs[1].Charge()*m_sW2);
  m_cR[1] = I*m_e/(2.*m_sW*m_cW)*(-2.*m_flavs[1].Charge()*m_sW2);
  InitAnalysis();
  msg->SetPrecision(16);
}

Z_Decay_Real::~Z_Decay_Real()
{
  FinishAnalysis(m_quark,m_lepton);
}

void Z_Decay_Real::InitAnalysis() {
  m_histos[string("shat")]	 = new Histogram(0,10.,10000.,10000);
  m_histos[string("m_omega")]	 = new Histogram(0,0,8000.*m_deltas,200);
  m_histos[string("Decay_mass_before")]	 = new Histogram(0,10.,500.,10000);
  m_histos[string("Decay_mass_after")]	 = new Histogram(0,10.,500.,10000);
  m_histos[string("mll_before")]	 = new Histogram(0,10.,500.,10000);
  m_histos[string("mll_after")]	 = new Histogram(0,10.,500.,10000);
  m_histos[string("Denominator_jac_L")]	 = new Histogram(0,10.,10000.,10000);
  m_histos[string("Denominator_jac_M")]	 = new Histogram(0,-10000.,10.,10000);
  m_histos[string("m_theta_closest")]	 = new Histogram(0,0.,10.,1000);
  m_histos[string("costheta")]	 = new Histogram(0,-1.,1.,200);
  m_histos[string("Eikonal")]	 = new Histogram(0,0.,1.E10,10000);
  m_histos[string("ME")]	 = new Histogram(0,0.,1.E10,10000);
  m_histos[string("ME_over_Eikonal")]	 = new Histogram(0,0.,1.E1,10000);
  // Enter res in this histogram, in 1st bin if w < born, in 2nd bin if w > born
  m_histos[string("Weights_above_below")]	 = new Histogram(0,0.,2.,2);
  m_histos[string("Numbers_above_below")]	 = new Histogram(0,0.,2.,2);
  // Overflow: If w > born, enter w-born into this histogram. If w < born, do not fill.
  m_histos[string("Overflow")]	 = new Histogram(0,0.,1.,1);
  m_histos[string("Overflow_numbers")]	 = new Histogram(0,0.,1.,1);
  m_histos[string("Full_weight")]	 = new Histogram(0,0.,1.E1,10000);
  m_histos[string("Born")]	 = new Histogram(0,0.,1.E4,10000);
  m_histos[string("E_gamma_gen")]    = new Histogram(0,0,10000.,1000);
  m_histos[string("E_gamma_acc")]    = new Histogram(0,0,10000.,1000);
  m_histos[string("W_J_L")]	 = new Histogram(0,0,10.,1000);
  m_histos[string("W_J_M")]	 = new Histogram(0,0,10.,1000);
  m_histos[string("m_u")]	 = new Histogram(0,-1.,2.,10000);
  m_histos[string("Thrown_away")]	 = new Histogram(0,0.,2.,2);
  m_histos[string("Cuts_Failed")]	 = new Histogram(0,0.,2.,2);
  m_histos[string("Mlgamma")]	 = new Histogram(0,0.,100.,1000);
  m_histos[string("E_gamma")]	 = new Histogram(0,0.,55.,550);
  m_histos[string("pT_gamma")]	 = new Histogram(0,0.,55.,550);
  m_histos[string("Mlgamma_small")]	 = new Histogram(0,0.,1.,1000);
  m_histos[string("E_gamma_small")]	 = new Histogram(0,0.,0.1,1000);
  m_histos[string("pT_gamma_small")]	 = new Histogram(0,0.,1.,1000);
}

void Z_Decay_Real::FinishAnalysis(const Flavour& fl, const Flavour& fl_lep) {
  Histogram * histo;
  string name;
  string flavname = fl.ShellName();
  string lepname = fl_lep.ShellName();
  for (map<string,Histogram *>::iterator 
     hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    histo->MPISync();
    if (m_deltas == 0.001) {
      name  = string("Debug_Analysis_0.001")+string("_")+lepname+string("/")+flavname+string("_")+hit->first+string(".dat");
    }
    else if (m_deltas == 0.01) {
      name  = string("Debug_Analysis_0.01")+string("_")+lepname+string("/")+flavname+string("_")+hit->first+string(".dat");
    }
    else if (m_deltas == 0.1) {
      name  = string("Debug_Analysis_0.1")+string("_")+lepname+string("/")+flavname+string("_")+hit->first+string(".dat");
    }
    else {
      name  = string("Debug_Analysis/")+flavname+string("_")+hit->first+string(".dat");
    }
    PRINT_VAR(name);
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}


double Z_Decay_Real::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  // None of these ever happen
  if ((p[2]+p[3]).Abs() < 50.) msg_Error() << "Born cut not respected. m_ll " << (p[2]+p[3]).Abs() << std::endl;
  if (p[2].PPerp() < 25. || p[3].PPerp() < 25.) msg_Error() << "Born cut not respected. p[2].pT = " << p[2].PPerp() << " , p[3].pT = " << p[3].PPerp() << std::endl;
  if (abs(p[2].Eta()) > 2.5 || p[3].Eta() > 2.5) msg_Error() << "Born cut not respected. p[2].eta = " << p[2].Eta() << " , p[3].eta = " << p[3].Eta() << std::endl;

  // Convert delta_s into cut on photon energy: Ey < delta_s*E_CMS/2 in partonic CMS
  m_omega = m_deltas*(p[0]+p[1]).Abs()/2.;
  // Set masses for decay setup
  m_masses[0] = (p[0]+p[1]).Abs();
  for (int i(1); i < 3; i++) {
    m_masses[i] = m_flavs[i].Mass();
  }
  // Determine maximum photon energy. To cover entire possible space, choose max such that u in [0,1]
  // u^2 = ((K0 - sqrt(M^2+K0^2))^2 - 4*ml^2)/(4*q^2) -> K0_max = M/2*(M/(2*ml)-2*ml/M)
  m_max = m_masses[0]/2.*(m_masses[0]/(2.*m_masses[1])-(2.*m_masses[1])/m_masses[0]);

  Complex res(0.,0.);
  // Boost particles into partonic CMS (also dipole CMS for lepton-antilepton dipole)
  Poincare CMS(p[0]+p[1]);
  Vec4D temp = p[0]+p[1];
  CMS.Boost(temp);
  m_moms0[0] = temp;
  temp = p[2];
  CMS.Boost(temp);
  m_moms0[1] = temp;
  temp = p[3];
  CMS.Boost(temp);
  m_moms0[2] = temp;
  // Rotate such that lepton lies in the z-direction, antilepton in -z-direction
  Poincare rot(m_moms0[1],Vec4D(0.,0.,0.,1.));
  rot.Rotate(m_moms0[1]);
  rot.Rotate(m_moms0[2]);






  bool passedcuts = true;
  double entry = 1.;
  double born = 0.;
  Setup_Corrected_Momenta(p);

  // Pseudorapidities and pTs have to be calculated in the lab frame
  // Boost into CMS of Z after mapping, rotate back such that lepton not on z-axis anymore
  // boost back into lab frame
  Poincare newCMS(m_moms[0]);
  Vec4D Zmom = m_moms[0];
  newCMS.Boost(Zmom);
  rot.RotateBack(Zmom);
  CMS.BoostBack(Zmom);
  Vec4D leptonmom = m_moms[1];
  newCMS.Boost(leptonmom);
  rot.RotateBack(leptonmom);
  CMS.BoostBack(leptonmom);  
  Vec4D antileptonmom = m_moms[2];
  newCMS.Boost(antileptonmom);
  rot.RotateBack(antileptonmom);
  CMS.BoostBack(antileptonmom);  
  Vec4D labphoton = m_moms[3];
  newCMS.Boost(labphoton);
  rot.RotateBack(labphoton);
  Vec4D labphoton_CMS = labphoton;
  Vec4D leptonmom_CMS = leptonmom;
  Vec4D antileptonmom_CMS = antileptonmom;
  CMS.BoostBack(labphoton);
  CMS.Boost(leptonmom_CMS);
  CMS.Boost(antileptonmom_CMS);
  Vec4D MomConservation = p[0]+p[1]-(leptonmom+antileptonmom+labphoton);
  // Generated photon energy is NOT in partonic CMS but in frame boosted wrt partonic CMS by k
  // hence need to boost to partonic CMS and check that it is above cut and below max energy photon
  // can carry away!
  m_histos[string("mll_before")]->Insert((m_moms0[1]+m_moms0[2]).Abs(),entry);  
  if (labphoton_CMS[0] < m_omega  || labphoton_CMS[0] > ((p[0]+p[1]).Abs2()-4.*sqr(m_masses[1]))/(2.*(p[0]+p[1]).Abs())
      ) {
    entry = 0.;
    res = 0.;
    passedcuts = false;
    m_histos[string("Thrown_away")]->Insert(0.5,1.);
  }
  // Check cuts after photon generation
  // Invariant masses can be checked in any frame
  // lepton invariant mass: > 50 GeV
  else if ((m_moms[1]+m_moms[2]).Abs() < 50.) {
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // lepton invariant mass: < sqrt(s)
  else if ((m_moms[1]+m_moms[2]).Abs() > (p[0]+p[1]).Abs()) {
    msg_Out() << "Huh? Lepton invariant mass " << (m_moms[1]+m_moms[2]).Abs() << " partonic COM energy " << (p[0]+p[1]).Abs() << "\n";
    // Invariant mass can be checked in any frame
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // lepton photon invariant mass: < sqrt(s)
  else if (((m_moms[1]+m_moms[2]+m_moms[3]).Abs()-(p[0]+p[1]).Abs()) > 1E-6*(p[0]+p[1]).Abs()) {
    msg_Out() << "Huh? Lepton-photon invariant mass " << (m_moms[1]+m_moms[2]+m_moms[3]).Abs() << " partonic COM energy " << (p[0]+p[1]).Abs() << "\n";
    res = 0.;
    entry = 0.;
    passedcuts = false;
  }
  // Check momentum conservation in lab frame
  else if (MomConservation[0] > 1.E-6*(p[0]+p[1]).Abs() ||
	   MomConservation[1] > 1.E-6*(p[0]+p[1]).Abs() ||
	   MomConservation[2] > 1.E-6*(p[0]+p[1]).Abs() ||
	   MomConservation[3] > 1.E-6*(p[0]+p[1]).Abs()) {
    msg_Out() << "Momentum not conserved after boost back. Diff " << MomConservation << "\n";
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }
  // Check lepton momenta are larger than 0
  else if (leptonmom[0] < 0. || antileptonmom[0] < 0.) {
    msg_Out() << "Lepton energy below 0.\nLepton " << leptonmom << "\nAntilepton "  << antileptonmom << std::endl;;
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }
  // Check photon momentum is larger than 0
  else if (labphoton[0] < 0.) {
    msg_Out() << "Photon energy below 0.\n " << labphoton << std::endl;;
    entry = 0.;
    res = 0.;
    passedcuts = false;
  }
  // Check other cuts: |eta(l)| < 2.5; pT(l) < 25 GeV
  else if (abs(leptonmom.Eta()) > 2.5 ||
	   abs(antileptonmom.Eta()) > 2.5 ||
	   leptonmom.PPerp() < 25. ||
	   antileptonmom.PPerp() < 25.) {
    res= 0.;
    entry = 0.;
    passedcuts = false;
  }

  // // Cuts on cos(theta) in CMS frame
  // else if ((Vec3D(labphoton_CMS)*Vec3D(antileptonmom))/(Vec3D(labphoton_CMS).Abs()*Vec3D(antileptonmom).Abs()) < -0.95 || (Vec3D(labphoton_CMS)*Vec3D(antileptonmom))/(Vec3D(labphoton_CMS).Abs()*Vec3D(antileptonmom).Abs()) > 0.95) {
  //   res = 0.;
  //   entry = 0.;
  //   passedcuts = false;
  // }
  
  if (passedcuts == false) {
    m_histos[string("Cuts_Failed")]->Insert(0.5,1.);
    m_histos["E_gamma"]->Insert(labphoton[0],0.);
    m_histos["pT_gamma"]->Insert(labphoton.PPerp(),0.);
    m_histos["Mlgamma"]->Insert((labphoton+antileptonmom).Abs(),0.);
    m_histos["E_gamma_small"]->Insert(labphoton[0],0.);
    m_histos["pT_gamma_small"]->Insert(labphoton.PPerp(),0.);
    m_histos["Mlgamma_small"]->Insert((labphoton+antileptonmom).Abs(),0.);
    m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(antileptonmom))/(Vec3D(labphoton).Abs()*Vec3D(antileptonmom).Abs()),0.);
  }

  if (passedcuts == true) {
    // Calculate born cross section - qedterm ~ |M^gamma|^2, Zterm ~ |M^Z|^2, intterm ~ interference
    m_histos[string("Cuts_Failed")]->Insert(1.5,1.);
    m_histos[string("Thrown_away")]->Insert(1.5,1.);
    double s(0.),t(0.);
    s=(p[0]+p[1]).Abs2();
    t=(p[0]-p[2]).Abs2();
    chi1  = m_kappa * s * (s-MZ2)/(sqr(s-MZ2) + GZ2*MZ2);
    chi2  = sqr(m_kappa * s)/(sqr(s-MZ2) + GZ2*MZ2);
    
    qedterm = (1+sqr(1.+2.*t/s)) * sqr(qf*qe);
    intterm = (1+sqr(1.+2.*t/s)) * 2.*(qf*qe*vf*ve) * chi1 + (1.+2.*t/s) * 4. * qe*qf*ae*af * chi1;
    Zterm = (1+sqr(1.+2.*t/s)) * (ae*ae+ve*ve) * (af*af+vf*vf) * chi2 + (1.+2.*t/s) * 8. * ae*ve*af*vf * chi2;
    born = qedterm+intterm+Zterm;
    // Writeout to compare MEs with Comix/Amegic
    // PRINT_VAR("\n  1  " << p[0][0] << "  " << p[0][1] <<  "  " << p[0][2] << "  " << p[0][3] << "\n  -1  " << p[1][0] << "  " << p[1][1] <<  "  " << p[1][2] << "  " << p[1][3] << "\n  11  " << p[2][0] << "  " << p[2][1] <<  "  " << p[2][2] << "  " << p[2][3] << "\n  -11  " << p[3][0] << "  " << p[3][1] <<  "  " << p[3][2] << "  " << p[3][3]);
    // PRINT_VAR(born*sqr(4.*M_PI*m_alpha*CouplingFactor(0,1))/3.);
    
    // Calculate weights - jac_L: weight due to LT into dipole CMS,
    // jac_M: weight due to momentum mapping, Uncorrected_Smod:
    // eikonal factor using uncorrected lepton momenta
    double weight_jac_L = Weight_Jac_L();
    double weight_jac_M = Weight_Jac_M();
    double Uncorrected_Smod = Smod_Uncorr();
    
    // Reweight each piece according to their respective decay ME
    // contribution - checked that these ar the same for all
    // contributions
    if (m_soft_approx == 0) {
      Complex qedres = GetBeta_1_1(0,0)/GetBeta_0_0(0,0)*qedterm*weight_jac_L*weight_jac_M/Uncorrected_Smod
	;
      Complex intres = (GetBeta_1_1(0,1)+GetBeta_1_1(1,0))/(GetBeta_0_0(0,1)+GetBeta_0_0(1,0))*intterm*weight_jac_L*weight_jac_M/Uncorrected_Smod
	;
      Complex Zres = GetBeta_1_1(1,1)/GetBeta_0_0(1,1)*Zterm*weight_jac_L*weight_jac_M/Uncorrected_Smod
	;
      // Safeguard against numerical issues
      if (Uncorrected_Smod < pow(10.,-25.) || m_angle_weight < pow(10.,-25.)) res = 0.;
      else res = m_alpha/(4.*M_PI*M_PI)*CouplingFactor(0,1)*log(m_max/m_omega)*m_angle_weight*(qedres+intres+Zres);
      if (std::isinf(res.real())) PRINT_VAR(log(m_max/m_omega) << " , " << m_angle_weight << " , " << GetBeta_1_1(1,1)/GetBeta_0_0(1,1) << " , " << weight_jac_L << " , " << weight_jac_M << " , " << Uncorrected_Smod);
    }
    else {
      res = Smod()/Uncorrected_Smod*weight_jac_L*weight_jac_M*(qedterm+intterm+Zterm);
      // res = (qedterm+intterm+Zterm);
    }
    m_histos["E_gamma"]->Insert(labphoton[0],res.real());
    m_histos["pT_gamma"]->Insert(labphoton.PPerp(),res.real());
    m_histos["Mlgamma"]->Insert((labphoton+antileptonmom).Abs(),res.real());
    m_histos["E_gamma_small"]->Insert(labphoton[0],res.real());
    m_histos["pT_gamma_small"]->Insert(labphoton.PPerp(),res.real());
    m_histos["Mlgamma_small"]->Insert((labphoton+antileptonmom).Abs(),res.real());
    m_histos["costheta"]->Insert((Vec3D(labphoton)*Vec3D(antileptonmom))/(Vec3D(labphoton).Abs()*Vec3D(antileptonmom).Abs()),res.real());

    if (MC_method == 0) {
      // Implement hit and miss - born is maximum value (i. e. when photon momentum is 0)
      double max = born;
      double randomthrow = ran->Get()*max;
      if (randomthrow > res.real()) { 
	entry = 0.;
	res = 0.;
      }
    }
    m_histos[string("mll_after")]->Insert((m_moms[1]+m_moms[2]).Abs(),entry);  
    double photon_lepton_angle = m_moms[3].Theta(m_moms[1]);
    double photon_antilepton_angle = m_moms[3].Theta(m_moms[2]);
    if (photon_lepton_angle < photon_antilepton_angle) {
      m_histos[string("m_theta_closest")]->Insert(photon_lepton_angle*91.1876/(2.*m_masses[1]),entry);
    }
    else {
      m_histos[string("m_theta_closest")]->Insert(photon_antilepton_angle*91.1876/(2.*m_masses[1]),entry);
    }
    // Fill histos
    m_histos[string("m_u")]->Insert(m_u,entry);
    m_histos[string("E_gamma_acc")]->Insert(labphoton[0],entry);
    m_histos[string("Full_weight")]->Insert((GetBeta_1_1(0,0)/GetBeta_0_0(0,0)).real()/Smod_Uncorr()*weight_jac_L*weight_jac_M,entry);
    m_histos[string("Born")]->Insert(born,entry);
    m_histos[string("Eikonal")]->Insert(Smod_Uncorr(),entry);  
    m_histos[string("ME")]->Insert((GetBeta_1_1(0,0)/GetBeta_0_0(0,0)).real(),entry);  
    m_histos[string("ME_over_Eikonal")]->Insert((GetBeta_1_1(0,0)/GetBeta_0_0(0,0)).real()/Smod_Uncorr(),entry);  
    m_histos[string("shat")]->Insert(sqrt((p[0]+p[1]).Abs2()),entry);  
    m_histos[string("Decay_mass_before")]->Insert(sqrt(m_moms0[0].Abs2()),entry);  
    m_histos[string("Decay_mass_after")]->Insert(sqrt(m_moms[0].Abs2()),entry);  
    m_histos[string("m_omega")]->Insert(m_omega,entry);  
    m_histos[string("W_J_L")]->Insert(weight_jac_L,entry);  
    m_histos[string("W_J_M")]->Insert(weight_jac_M,entry);  
    if (real(res/born) > 1.) {
      m_histos[string("Overflow")]->Insert(0.5,sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*res.real());
      m_histos[string("Weights_above_below")]->Insert(1.5,sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*res.real());
      m_histos[string("Overflow_numbers")]->Insert(0.5,entry);
      m_histos[string("Numbers_above_below")]->Insert(1.5,entry);
    }
    else {
      m_histos[string("Weights_above_below")]->Insert(0.5,sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*res.real());
      m_histos[string("Numbers_above_below")]->Insert(0.5,entry);
    }
    // Check what causes Nan if they appear
    if (IsNan(res.real())) {
      PRINT_VAR("\n" << m_moms[0] << "\n" << m_moms[1] << "\n" << m_moms[2] << "\n" << m_moms[3]);
      PRINT_VAR(res << " , " << weight_jac_L << " , " << weight_jac_M << " , " << Smod_Uncorr() << " , " << (GetBeta_1_1(0,0)/GetBeta_0_0(0,0)));
    }
  }
  // 1/3 factor due to quark colour averaging
  if (MC_method == 0) {
    // Hit and miss
    return sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*entry*born;
  }
  else if (MC_method == 1) {
    // Entering probability directly
    return sqr(4.*M_PI*m_alpha)*CouplingFactor(0,2)*1./3.*res.real();
  }
}

void Z_Decay_Real::Setup_Corrected_Momenta(const ATOOLS::Vec4D_Vector& p) {
  // Generate photon
  bool reject = true;
  double m = 0.5*(m_masses[1]+m_masses[2]);
  double m2 = sqr(m);
  double M2 = sqr(m_masses[0]);
  while (reject) {
    // Generate photon
    m_moms[3] = Generate_One_Photon(m_moms0[1],m_moms0[2]);
    
    // determine u in mapping procedure
    double en = sqrt(M2+Vec3D(m_moms[3])*Vec3D(m_moms[3]));
    double numerator = sqr(m_moms[3][0]-en) - 4.*m2;
    Vec3D q = Vec3D(m_moms0[1]);
    double denominator = 4.*q*q;
    // Veto against very small and negative u to prevent numerical issues
    if (numerator/denominator > 0.) reject = false;
    m_u = sqrt(numerator/denominator);
  }

  // Rescale lepton momenta (taking into account masses correctly)
  Vec3D mom3 = m_u*Vec3D(m_moms0[1]);
  double momE = sqrt(m2+mom3.Sqr());
  m_moms[1] = Vec4D(momE,mom3);
  mom3 = m_u*Vec3D(m_moms0[2]);
  momE = sqrt(m2+mom3.Sqr());
  m_moms[2] = Vec4D(momE,mom3);
  // Compensate momenta with Z-momentum
  m_moms[0] = m_moms[1]+m_moms[2]+m_moms[3];
  // make sure Z is still on-shell
  if (abs(m_moms[0].Abs()-m_moms0[0].Abs()) > 1.E-6*m_moms[0][0]) {
    msg_Out() << METHOD << "decay mass moved off-shell! before " << m_moms0[0].Abs() << " after " << m_moms[0].Abs() << " diff " << abs(m_moms0[0].Abs()-m_moms[0].Abs()) << " \n";
    msg_Out() << METHOD << "Z Energy before " << m_moms0[0][0] << " after " << m_moms[0][0] << " photon momentum " << m_moms[3] << " u " << m_u << "\n";
  }
}

double Z_Decay_Real::Weight_Jac_L() {
  // Energy of charged particles before mapping
  double QC0 = (m_moms0[1]+m_moms0[2])[0];
  // Energy of charged particles after mapping
  double PC0 = (m_moms[1]+m_moms[2])[0];
  // In dipole rest frame hence E = M
  double mMQ = QC0;
  double mMP = PC0;
  // No additional neutral particles in process
  double QN0 = 0.;
  double PN0 = 0.;
  double K0 = m_moms[3][0];
  if (K0 < m_deltas*(m_moms0[1]+m_moms0[2]).Abs()/2.) {
    msg_Out() << METHOD << ": This should not happen. Photon energy less than cut allows.\n";
    THROW(fatal_error, "cannot continue");
  }
  m_histos[string("Denominator_jac_L")]->Insert(PC0+PN0+K0);  
  
  return pow(mMP,3.)/pow(mMQ,3.)*(QC0+QN0)/(PC0+PN0+K0);
}

double Z_Decay_Real::Weight_Jac_M() {
  double sumq     = 0.;
  double sumpq     = 0.;
  double prod     = 1.;
  double N = 2;
  for (int i(1); i < 3; i++) {
    sumq = sumq - Vec3D(m_moms0[i]).Sqr()/m_moms0[i][0];
    sumpq = sumpq - Vec3D(m_moms0[i])*Vec3D(m_moms[i])/m_moms[i][0];
    prod = prod*m_moms0[i][0]/m_moms[i][0];
  }
  m_histos[string("Denominator_jac_M")]->Insert(sumpq);  
  return pow(m_u,3.*N-4.) * sumq/sumpq * prod;
}

Complex Z_Decay_Real::InfraredSubtractedME_0_0(const int& b) {
  // Calculate leading order helicity amplitude
  // b = 0: photon mediates, b = 1: Z mediates
  Vec4C epsZ = Polarization_Vector(m_moms0[0])[m_spins0[0]];
  XYZFunc XYZ(3,m_moms0,m_flavs,false);
  return  XYZ.X(1,m_spins0[1],epsZ,2,m_spins0[2],m_cR[b]*sqrt(CouplingFactor(0,1)),m_cL[b]*sqrt(CouplingFactor(0,1)));
}

double Z_Decay_Real::Generate_Energy(const double& omega, const double& max) {
  // Generate energy following 1/E distribution
  return omega*pow(max/omega,ran->Get());
}

void Z_Decay_Real::Generate_Photon_Angle(const Vec4D& p1, const Vec4D& p2) {
  // No need to boost/rotate here as already in right frame

  // Calculate beta-factors
  double b1 = Vec3D(p1).Abs()/p1[0];
  double b2 = Vec3D(p2).Abs()/p2[0];
  // Generate dipole angles
  GenerateDipoleAngle(b1,b2);
  // Generate null vector from angles
  GenerateNullVector();
}

void Z_Decay_Real::GenerateDipoleAngle(const double& b1, const double& b2) {
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
  m_angle_weight = 2.*M_PI*weight*((1.+b1*b2)/(b1+b2)*log((1.+b1)*(1.+b2)/((1.-b1)*(1.-b2)))-2.);

  m_theta = acos(m_c);
  if (m_c < -1. || m_c > 1.) {
    msg_Out() << "Encountered c outside [-1,1], c = " << m_c << ". Theta = " << m_theta << "\n";
  }
  m_phi   = 2.*M_PI*ran->Get();
}

void Z_Decay_Real::GenerateNullVector() {
  m_dir = Vec4D(1., sin(m_theta)*cos(m_phi), sin(m_theta)*sin(m_phi), cos(m_theta));
}

Vec4D Z_Decay_Real::Generate_One_Photon(const Vec4D& p1, const Vec4D& p2) {
  // Generate photon energy and angle and return photon momentum
  double E = Generate_Energy(m_omega,m_max);
  if (E < m_deltas*(p1+p2).Abs()/2.) {
    msg_Out() << METHOD << ": This should not happen. Photon energy less than cut allows.\n";
    THROW(fatal_error, "cannot continue");
  }
  m_histos[string("E_gamma_gen")]->Insert(E);

  Generate_Photon_Angle(p1,p2);
  return E*m_dir;
}

Complex Z_Decay_Real::GetBeta_0_0(const int& b1, const int& b2) {
  // Sum up helicity amplitudes for born decay cross section
  // b1, b2 determine which particle mediates - 0: photon, 1: Z
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
        m_spins0[0] = k;
        m_spins0[1] = j;
        m_spins0[2] = i;
        Complex M_0_0_b1 = InfraredSubtractedME_0_0(b1);
        Complex M_0_0_b2 = InfraredSubtractedME_0_0(b2);
        sum = sum + (M_0_0_b1*conj(M_0_0_b2))// .real()
	  ;
      }
    }
  }
  // spin average over initial state
  sum = (1./3.)*sum;
  return sum;
}

double Z_Decay_Real::Smod() {
  // Smod calculated in case of one additional photon
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double Z_Decay_Real::Smod_Uncorr() {
  // Smod calculated in case of one additional photon using uncorrected momenta
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms0[1];
  Vec4D pj  = m_moms0[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double Z_Decay_Real::Smod_Norm() {
  // Smod calculated in case of one additional photon using uncorrected momenta
  Vec4D k   = m_moms[3]/m_moms[3][0];
  Vec4D pi  = m_moms0[1];
  Vec4D pj  = m_moms0[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha*CouplingFactor(0,1)/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex Z_Decay_Real::InfraredSubtractedME_1_05(const int& b) {
  // Calculate real ME
  // Copy momenta just to be sure
  Vec4D real_moms[8];
  for (size_t j = 0; j < 4; j++) real_moms[j] = m_moms[j];
  Vec4C epsZ   = Polarization_Vector(real_moms[0])[m_spins[0]];
  Vec4C epsP   = conj(Polarization_Vector(real_moms[3])[m_spins[3]]);
  Vec4D pa     = real_moms[1]+real_moms[3];	    // fermion propagator momenta
  Vec4D pb     = real_moms[2]+real_moms[3];
  double m     = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
  real_moms[4]    = real_moms[5] = pa;	       // enter those into xi*real_moms
  real_moms[6]    = real_moms[7] = pb;
  m_flavs[4]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[7] = m_flavs[2];
  // Note the massmode - internal fermions have mass = \sqrt(p^2)
  XYZFunc XYZ(8,real_moms,m_flavs,false);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  // emission off fermion (a) and (b)
  for (unsigned int s=0; s<=1; s++) {
    r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)
          *XYZ.X(4,s,epsZ,2,m_spins[2],m_cR[b],m_cL[b]);
    r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)
          *XYZ.X(5,s,epsZ,2,m_spins[2],m_cR[b],m_cL[b]);
    r3 += XYZ.X(1,m_spins[1],epsZ,6,s,m_cR[b],m_cL[b])
          *XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    r4 += XYZ.X(1,m_spins[1],epsZ,7,s,m_cR[b],m_cL[b])
          *XYZ.X(7,s,epsP,2,m_spins[2],1.,1.);
  }
  // add prefactors
  r1 *= m_e*CouplingFactor(0,1)/(2.*(pa*pa-m*m))*(1.+m/sqrt(pa*pa));
  r2 *= m_e*CouplingFactor(0,1)/(2.*(pa*pa-m*m))*(1.-m/sqrt(pa*pa));
  r3 *= -m_e*CouplingFactor(0,1)/(2.*(pb*pb-m*m))*(1.-m/sqrt(pb*pb));
  r4 *= -m_e*CouplingFactor(0,1)/(2.*(pb*pb-m*m))*(1.+m/sqrt(pb*pb));
  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = m_flavs[6] = m_flavs[7] = Flavour(kf_none);
  return (r1+r2+r3+r4);
}

Complex Z_Decay_Real::GetBeta_1_1(const int& b1, const int& b2) {
  // Assemble real matrix element squared
  Complex sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
        for (unsigned int l=0; l<=1; l++) {     // spin gamma
          m_spins[0] = k;
          m_spins[1] = j;
          m_spins[2] = i;
          m_spins[3] = l;
	  Complex M_1_05_b1 = InfraredSubtractedME_1_05(b1);
	  Complex M_1_05_b2 = InfraredSubtractedME_1_05(b2);
          sum = sum + (M_1_05_b1*conj(M_1_05_b2))// .real()
	    ;
        }
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  sum = 1./(16.*M_PI*M_PI*M_PI)*sum;
  return sum;
}

DECLARE_TREEME2_GETTER(Z_Decay_Real,
		       "Z_Decay_Real")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,Z_Decay_Real>::
operator()(const Process_Info &pi) const
{
  Default_Reader reader;
  if (reader.GetValue<int>("EXTRAXS_Z_Decay_Real",0) != 1) return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
      fl[2].IsLepton()  && fl[2]==fl[3].Bar()) {
    if (pi.m_maxcpl[1]==3 && pi.m_mincpl[1]==3) {
      return new Z_Decay_Real(pi,fl);
    }
  }
  return NULL;
}


