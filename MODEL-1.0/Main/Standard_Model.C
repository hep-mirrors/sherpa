#include "Standard_Model.H"
#include "Running_AlphaQED.H"
#include "Running_AlphaS.H"
#include "Running_Fermion_Mass.H"
#include "Effective_Higgs_Coupling.H"
#include "Message.H"
#include "Hdecay_Fortran_Interface.H"

using namespace MODEL;
using namespace ATOOLS;


Standard_Model::Standard_Model(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the Standard Model from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();

  ReadInFile();
}

void Standard_Model::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  FixEWParameters();  
  FixCKM();
  p_constants->insert(std::make_pair(std::string("Yukawa_e"), 
				     p_dataread->GetValue<double>("YUKAWA_E",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_mu"), 
				     p_dataread->GetValue<double>("YUKAWA_MU",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_tau"), 
				     p_dataread->GetValue<double>("YUKAWA_TAU",Flavour(kf::tau).PSMass())));
  p_constants->insert(std::make_pair(std::string("Yukawa_d"), 
				     p_dataread->GetValue<double>("YUKAWA_D",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_u"), 
				     p_dataread->GetValue<double>("YUKAWA_U",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_s"), 
				     p_dataread->GetValue<double>("YUKAWA_S",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_c"), 
				     p_dataread->GetValue<double>("YUKAWA_C",0.)));
  p_constants->insert(std::make_pair(std::string("Yukawa_b"), 
				     p_dataread->GetValue<double>("YUKAWA_B",Flavour(kf::b).PSMass())));
  p_constants->insert(std::make_pair(std::string("Yukawa_t"), 
				     p_dataread->GetValue<double>("YUKAWA_T",Flavour(kf::t).PSMass())));

  int    order_alphaS	= p_dataread->GetValue<int>("ORDER_ALPHAS",0);
  double alphaS         = p_dataread->GetValue<double>("ALPHAS(MZ)",0.1188);
  double alphaS_default = p_dataread->GetValue<double>("ALPHAS(default)",alphaS);
  double MZ2            = sqr((*p_constants)[std::string("MZ")]);

  as = new Running_AlphaS(alphaS,MZ2,order_alphaS);
  as->SetDefault(alphaS_default);

  p_constants->insert(std::make_pair(std::string("alpha_S(MZ)"),alphaS));
  p_functions->insert(std::make_pair(std::string("alpha_S"),as));

  Running_Fermion_Mass * md   = new Running_Fermion_Mass(Flavour(kf::d),
							 ScalarConstant(std::string("Yukawa_d")),as);
  Running_Fermion_Mass * mu   = new Running_Fermion_Mass(Flavour(kf::u),
							 ScalarConstant(std::string("Yukawa_u")),as);
  Running_Fermion_Mass * ms   = new Running_Fermion_Mass(Flavour(kf::s),
							 ScalarConstant(std::string("Yukawa_s")),as);
  Running_Fermion_Mass * mc   = new Running_Fermion_Mass(Flavour(kf::c),
							 ScalarConstant(std::string("Yukawa_c")),as);
  Running_Fermion_Mass * mb   = new Running_Fermion_Mass(Flavour(kf::b),
							 ScalarConstant(std::string("Yukawa_b")),as);
  Running_Fermion_Mass * mt   = new Running_Fermion_Mass(Flavour(kf::t),
							 ScalarConstant(std::string("Yukawa_t")),as);
  Running_Fermion_Mass * me   = new Running_Fermion_Mass(Flavour(kf::e),
							 ScalarConstant(std::string("Yukawa_e")),as);
  Running_Fermion_Mass * mmu  = new Running_Fermion_Mass(Flavour(kf::mu),
							 ScalarConstant(std::string("Yukawa_mu")),as);
  Running_Fermion_Mass * mtau = new Running_Fermion_Mass(Flavour(kf::tau),
							 ScalarConstant(std::string("Yukawa_tau")),as);
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::d).Name()),md));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::u).Name()),mu));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::s).Name()),ms));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::c).Name()),mc));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::b).Name()),mb));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::t).Name()),mt));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::e).Name()),me));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::mu).Name()),mmu));
  p_functions->insert(std::make_pair(std::string("m")+std::string(Flavour(kf::tau).Name()),mtau));

  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 
  double eh=2./3.;
  if (p_dataread->GetValue<int>("FINITE_TOP_MASS",0)==1) {
    double hm=Flavour(kf::h).Mass();
    Effective_Higgs_Coupling ehc(hm);
    eh = ehc.GetFermionContribution(Flavour(kf::t).Mass());
  }
  p_constants->insert(std::make_pair(std::string("Higgs_gg_fac"),eh));
}


void Standard_Model::FixEWParameters() {
  double MW,MZ,MH,alphaQED,sin2thetaW,cos2thetaW,vev,lambdaH,GF;
  m_ewscheme = p_dataread->GetValue<int>("EW_SCHEME",0);
  switch (m_ewscheme) {
  case 1:
    // SM parameters given by alphaQED, M_W, M_Z, M_H
    alphaQED   = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    MW         = Flavour(kf::W).Mass();
    MZ         = Flavour(kf::Z).Mass();
    MH         = Flavour(kf::h).Mass();
    cos2thetaW = sqr(MW/MZ);
    sin2thetaW = 1.-cos2thetaW;
    vev        = 2.*MW*sqrt(sin2thetaW/(4.*M_PI*alphaQED));
    lambdaH    = 2.*sqr(MH/vev);
    break;
  case 2:
    // SM parameters given by alphaQED, sinthetaW, v, M_H
    alphaQED   = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    vev        = p_dataread->GetValue<double>("VEV",246.);
    sin2thetaW = p_dataread->GetValue<double>("SIN2THETAW",0.23);
    cos2thetaW = 1.-sin2thetaW;
    MW         = vev/2.*sqrt((4.*M_PI*alphaQED)/sin2thetaW);
    MZ         = vev/2.*sqrt((4.*M_PI*alphaQED)*(1/sin2thetaW+1/cos2thetaW));
    MH         = p_dataread->GetValue<double>("MH",120.);
    lambdaH    = 2.*sqr(MH/vev);
  default:
    // all SM parameters given explicitly
    alphaQED   = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    MW         = Flavour(kf::W).Mass();
    MZ         = Flavour(kf::Z).Mass();
    MH         = Flavour(kf::h).Mass();
    sin2thetaW = p_dataread->GetValue<double>("SIN2THETAW",0.23);
    cos2thetaW = 1.-sin2thetaW;
    vev        = p_dataread->GetValue<double>("VEV",246.);
    lambdaH    = p_dataread->GetValue<double>("LAMBDA",0.47591);
    break;
  }
  aqed                    = new Running_AlphaQED(alphaQED,sqr(MZ));
  double alphaQED_default = 1./p_dataread->GetValue<double>("1/ALPHAQED(default)",1./(*aqed)(sqr(MZ)));
  aqed->SetDefault(alphaQED_default);

  p_functions->insert(std::make_pair(std::string("alpha_QED"),aqed));

  GF = sqrt(2.)*(*aqed)(sqr(Flavour(kf::mu).PSMass()))*M_PI/(2.*sin2thetaW*sqr(MW));

  p_constants->insert(std::make_pair(std::string("alpha_QED(0)"),alphaQED));
  p_constants->insert(std::make_pair(std::string("sin2_thetaW"), sin2thetaW));
  p_constants->insert(std::make_pair(std::string("cos2_thetaW"), cos2thetaW));
  p_constants->insert(std::make_pair(std::string("vev"),         vev));
  p_constants->insert(std::make_pair(std::string("MW"),          MW));
  p_constants->insert(std::make_pair(std::string("MZ"),          MZ));
  p_constants->insert(std::make_pair(std::string("GammaW"),      Flavour(kf::W).Width()));
  p_constants->insert(std::make_pair(std::string("GammaZ"),      Flavour(kf::Z).Width()));
  p_constants->insert(std::make_pair(std::string("MH"),          MH));
  p_constants->insert(std::make_pair(std::string("lambdaH"),     lambdaH));
  p_constants->insert(std::make_pair(std::string("GF"),          GF));
}

void Standard_Model::FixCKM() {
  CMatrix CKM(3);

  for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) CKM[i][j] = CKM[j][i] = Complex(0.,0.);
    CKM[i][i] = Complex(1.,0.);
  }
  
  double Cabibbo,A,rho,eta;
  m_ckmorder     = p_dataread->GetValue<int>("CKMORDER",0);  
  if (m_ckmorder>0) {
    Cabibbo    = p_dataread->GetValue<double>("CABIBBO",0.22);
    CKM[0][0] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[1][1] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[0][1] += Cabibbo * Complex( 1.,0.);
    CKM[1][0] += Cabibbo * Complex(-1.,0.);
  }
  if (m_ckmorder>1) {
    A          = p_dataread->GetValue<double>("A",0.8);
    CKM[1][2] += A*sqr(Cabibbo)  * Complex( 1.,0.);
    CKM[2][1] += A*sqr(Cabibbo)  * Complex(-1.,0.);
  }
  if (m_ckmorder>2) {
    eta        = p_dataread->GetValue<double>("ETA",0.5);
    rho        = p_dataread->GetValue<double>("RHO",0.5);
    CKM[1][2] += A*sqr(Cabibbo) * Complex(rho,-eta);
    CKM[2][1] += A*sqr(Cabibbo) * Complex(1.-rho,-eta);
  }
  p_matrices->insert(std::make_pair(std::string("CKM"),CKM));
}


bool Standard_Model::RunSpectrumGenerator() {
  m_spectrum = p_dataread->GetValue<int>("GENERATOR_ON",0);
  if (m_spectrum) {
    m_generator = p_dataread->GetValue<std::string>("HIGGS_GENERATOR",std::string("Hdecay"));
    if (m_generator==std::string("Hdecay")) {
      p_spectrumgenerator = new HDECAY::Hdecay_Fortran_Interface(p_dataread,this);
      p_spectrumgenerator->Run(std::string("SM"));
      return 1;
    }
    
    msg.Error()<<"Error in Standard_Model::RunSpectrumGenerator."<<std::endl
	       <<"   Unknown spectrum generator : "<<m_generator<<" use internal solution."<<std::endl;
    return 0;
  }
  return 1;
}

bool Standard_Model::FillDecay(ATOOLS::Decay_Table * dt) 
{
  if (m_generator==std::string("Hdecay")) {
    return p_spectrumgenerator->FillDecay(dt);
  }
  return 0;
}

