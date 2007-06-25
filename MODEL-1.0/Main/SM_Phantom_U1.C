#include "SM_Phantom_U1.H"
#include "Standard_Model.H"
#include "Message.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;


SM_Phantom_U1::SM_Phantom_U1(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the SM_Phantom_U1 from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM_Phantom_U1");

  Model_Base * SM = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(SM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(SM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(SM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(SM->GetComplexMatrices()));

  ReadInFile();
  if (!SanityChecks()) {
    msg.Error()<<"Potential Error in "<<METHOD<<":"<<endl
	       <<"   Sanity checks not passed."<<endl
	       <<"   Continue and hope for the best."<<endl;
  }
  ConstructTriangleFactors();
  FillMasses();
  FillWidths();
}

SM_Phantom_U1::~SM_Phantom_U1() 
{
  if (!m_dectables.empty()) {
    for (map<Flavour, Decay_Table *>::iterator pdit=m_dectables.begin();
	 pdit!=m_dectables.end();pdit++) {
      delete pdit->second; pdit->second=NULL; 
    }
    m_dectables.clear();
  }
}

void SM_Phantom_U1::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  p_constants->insert(make_pair(string("Tan(Beta)"),    
				p_dataread->GetValue<double>("Tan(Beta)",1.)));
  p_constants->insert(make_pair(string("Tan(Theta)"),    
				p_dataread->GetValue<double>("Tan(Theta)",0.)));
  p_constants->insert(make_pair(string("M_H1"),    
				p_dataread->GetValue<double>("M_H1",-1.)));
  p_constants->insert(make_pair(string("M_H2"),    
				p_dataread->GetValue<double>("M_H2",-1.)));
  p_constants->insert(make_pair(string("M_Z'"),    
				p_dataread->GetValue<double>("M_Z'",-1.)));
  p_constants->insert(make_pair(string("g'_1"),    
				p_dataread->GetValue<double>("g'_1",0.)));

  CMatrix HiggsMix(2);
  HiggsMix[0][0] = HiggsMix[1][1] = sqrt(1./(1.+sqr(ScalarConstant(string("Tan(Theta)")))));
  HiggsMix[1][0] = sqrt(1.-sqr(abs(HiggsMix[1][1])));
  HiggsMix[0][1] = -HiggsMix[1][0];


  p_matrices->insert(std::make_pair(std::string("HiggsMix"),HiggsMix));
}

void SM_Phantom_U1::FillMasses() {
  Flavour flav;
  flav = Flavour(kf::h);
  flav.SetOn(false);
  // Switch off standard model Higgs sector completely

  flav = Flavour(kf::h0);
  flav.SetMass(ScalarConstant(string("M_H1")));
  flav.SetMassOn(true);
  flav = Flavour(kf::H0);
  flav.SetMass(ScalarConstant(string("M_H2")));
  flav.SetMassOn(true);
  flav = Flavour(kf::A0);
  flav.SetMass(0.);
  flav.SetMassOn(true);
  flav = Flavour(kf::ZPrime);
  flav.SetMass(ScalarConstant(string("M_Z'")));
  flav.SetMassOn(true);
}


void SM_Phantom_U1::FillWidths() {  
  Decay_Table   * dectable;
  Decay_Channel * channel;
  double width, mix;

  set<Flavour>    flouts;
  FillPotentialDecayProducts(flouts);

  Flavour flin(kf::h0);
  dectable = new Decay_Table(flin);
  m_dectables[flin] = dectable;
  dectable->SetWidthGenerator(string("MODEL::SM_Phantom (Sherpa)"));
  mix = sqr(abs(ComplexMatrixElement(string("HiggsMix"),0,0)));
  for (set<Flavour>::iterator flout=flouts.begin();flout!=flouts.end();flout++) {
    if ((*flout)==Flavour(kf::photon))     width = H2PhotonDecay(flin);
    else if ((*flout)==Flavour(kf::gluon)) width = H2GluonDecay(flin);
    else if ((*flout).IsFermion())         width = H2FDecay(flin,(*flout));
    else if ((*flout).IsScalar())          width = H2SDecay(flin,(*flout));
    else if ((*flout).IsVector())          width = H2VDecay(flin,(*flout));
    else                                   width = 0.;
    if (width>0.) {
      channel = new Decay_Channel(flin);
      channel->AddDecayProduct((*flout));
      channel->AddDecayProduct((*flout).Bar());
      channel->SetWidth(width*mix);
      dectable->AddDecayChannel(channel);
    }
  }
  flin.SetWidth(dectable->TotalWidth());
  if (msg.LevelIsInfo() || msg.LevelIsTracking() || msg.LevelIsDebugging()) {
    dectable->Output();
  }


  flin = Flavour(kf::H0);
  dectable = new Decay_Table(flin);
  m_dectables[flin] = dectable;
  dectable->SetWidthGenerator(string("MODEL::SM_Phantom (Sherpa)"));
  mix = sqr(abs(ComplexMatrixElement(string("HiggsMix"),1,0)));
  for (set<Flavour>::iterator flout=flouts.begin();flout!=flouts.end();flout++) {
    if ((*flout)==Flavour(kf::photon))     width = H2PhotonDecay(flin);
    else if ((*flout)==Flavour(kf::gluon)) width = H2GluonDecay(flin);
    else if ((*flout).IsFermion())         width = H2FDecay(flin,(*flout));
    else if ((*flout).IsScalar())          width = H2SDecay(flin,(*flout));
    else if ((*flout).IsVector())          width = H2VDecay(flin,(*flout));
    else                                   width = 0.;
    if (width>0.) {
      channel = new Decay_Channel(flin);
      channel->AddDecayProduct((*flout));
      channel->AddDecayProduct((*flout).Bar());
      channel->SetWidth(width*mix);
      dectable->AddDecayChannel(channel);
    }
  }
  flin.SetWidth(dectable->TotalWidth());
  if (msg.LevelIsInfo() || msg.LevelIsTracking() || msg.LevelIsDebugging()) {
    dectable->Output();
  }
}

void SM_Phantom_U1::FillPotentialDecayProducts(set<Flavour> & flouts) {
  flouts.insert(Flavour(kf::photon));
  flouts.insert(Flavour(kf::gluon));
  flouts.insert(Flavour(kf::tau));
  flouts.insert(Flavour(kf::b));
  flouts.insert(Flavour(kf::t));
  flouts.insert(Flavour(kf::W));
  flouts.insert(Flavour(kf::Z));
  flouts.insert(Flavour(kf::h0));
  flouts.insert(Flavour(kf::H0));
  flouts.insert(Flavour(kf::A0));
}

double SM_Phantom_U1::H2FDecay(const Flavour & flin,const Flavour & flout) {
  if (flin.Mass()<2.*flout.Mass()) return 0.;
  double pref  = 1./(8.*M_PI);
  double beta  = sqrt(1.-4.*sqr(flout.Mass()/flin.Mass()));
  double massF = ScalarFunction((string("m")+flout.Name()),sqr(flin.Mass()));
  if (massF==0.) massF = flout.Mass();
  double width = pref*pow(beta,3.)*sqr(massF/ScalarConstant(string("vev")))*flin.Mass();
  if (flout.Strong()) width*=3.;
  //  width *= 3.*(1.+5.67/M_PI*ScalarFunction(string("alpha_S"),sqr(flin.Mass())));
  return width;
}

double SM_Phantom_U1::H2VDecay(const Flavour & flin,const Flavour & flout) {
  double pref, x=sqr(flout.Mass()/flin.Mass()), R, width;
  if (flin.Mass()>2.*flout.Mass()) {
    pref  = ScalarConstant(string("GF"))/(16.*sqrt(2.)*M_PI);
    R     = (1.-4.*x+12.*x*x)*sqrt(1.-4.*x);
    width = pref*R*pow(flin.Mass(),3.);
    if (flout==Flavour(kf::W)) width *= 2.;
  }
  else {
    pref  = 3.*sqr(ScalarConstant(string("GF")))/(16.*pow(M_PI,3.));
    R     = 
      3.*(1.-8.*x+20.*x*x)/sqrt(4.*x-1.)*acos((3.*x-1)/(2.*pow(x,1.5))) -
      (1.-x)/(2.*x)*(2.-13.*x+47.*x*x) - 3.*(1.-6.*x+4.*x*x)/2.*log(x);
    if (flout==Flavour(kf::Z)) 
      R  *=
	7./12. - 10./9.*ScalarConstant(string("sin2_thetaW"))+
	40.*sqr(ScalarConstant(string("sin2_thetaW")))/27.;
    width = pref*R*pow(flout.Mass(),4.)*flin.Mass();
  }
  return width;
}

double SM_Phantom_U1::H2SDecay(const Flavour & flin,const Flavour & flout) {
  if (flin.Mass()<2.*flout.Mass()) return 0.;
  double pref(0.),beta = sqrt(1.-4.*sqr(flout.Mass()/flin.Mass())), mix;
  if (flout==Flavour(kf::A0)) {
    pref = sqr(ScalarConstant(string("Tan(Beta)"))/ScalarConstant(string("vev")))*
      pow(flin.Mass(),5.);
    if (flin==Flavour(kf::h0)) {
      pref *= sqr(abs(ComplexMatrixElement(string("HiggsMix"),0,1))/
		  abs(ComplexMatrixElement(string("HiggsMix"),0,0)));
    }
    else if (flin==Flavour(kf::H0)) {
      pref *= sqr(abs(ComplexMatrixElement(string("HiggsMix"),1,1))/
		  abs(ComplexMatrixElement(string("HiggsMix"),1,0)));
    }
  }
  else if (flin==Flavour(kf::h0) && flout==Flavour(kf::H0)) {
    mix  = sqr(abs(ComplexMatrixElement(string("HiggsMix"),1,0)*
		   (ComplexMatrixElement(string("HiggsMix"),0,1)+
		    ComplexMatrixElement(string("HiggsMix"),1,1)*
		    ScalarConstant(string("Tan(Beta)")))));
    pref = sqr((sqr(flin.Mass())+2.*sqr(flout.Mass()))/ScalarConstant(string("vev")));
    pref *= mix*flin.Mass();
  }
  else if (flin==Flavour(kf::H0) && flout==Flavour(kf::h0)) {
    mix  = sqr(abs(ComplexMatrixElement(string("HiggsMix"),0,0)*
		   (ComplexMatrixElement(string("HiggsMix"),0,0)+
		    ComplexMatrixElement(string("HiggsMix"),1,0)*
		    ScalarConstant(string("Tan(Beta)")))));
    pref = sqr((sqr(flin.Mass())+2.*sqr(flout.Mass()))/ScalarConstant(string("vev")));
    pref *= mix*flin.Mass();
  }
  return pref*beta/(2.*16.*M_PI*sqr(flin.Mass()));
}

double SM_Phantom_U1::H2GluonDecay(const Flavour & flin) {
  double pref  = sqr(1./ScalarConstant(std::string("vev")))/(16.*2.*pow(M_PI,3.));
  if (flin==Flavour(kf::h0)) 
    pref *= sqr(ScalarConstant(std::string("Higgs_GG_eff_h")));
  else if (flin==Flavour(kf::H0)) 
    pref *= sqr(ScalarConstant(std::string("Higgs_GG_eff_H")));
  else
    pref *= 0.;
  double alp2  = sqr(ScalarFunction(string("alpha_S"),sqr(flin.Mass())));
  double width = pref*alp2*pow(flin.Mass(),3.);
  return width;
}

double SM_Phantom_U1::H2PhotonDecay(const Flavour & flin) {
  return 0.;
}

bool SM_Phantom_U1::SanityChecks() {
  if (ScalarConstant(string("Tan(Beta)"))<1.e-6 ||
      ScalarConstant(string("M_H1"))<0.       || 
      ScalarConstant(string("M_H2"))<0.) return false;
  return true;
}

void SM_Phantom_U1::ConstructTriangleFactors() {
  if (p_dataread->GetValue<int>("FINITE_TOP_MASS",0)==1) {
    double massratio=Flavour(kf::h0).Mass()/Flavour(kf::t).Mass();
    p_constants->insert(std::make_pair(std::string("Higgs_GG_eff_h"),
				       GetFermionContribution(massratio)));
    massratio = Flavour(kf::H0).Mass()/Flavour(kf::t).Mass();
    p_constants->insert(std::make_pair(std::string("Higgs_GG_eff_H"),
				       GetFermionContribution(massratio)));
  }
  else {
    p_constants->insert(std::make_pair(std::string("Higgs_GG_eff_h"),2./3.));
    p_constants->insert(std::make_pair(std::string("Higgs_GG_eff_H"),2./3.));
  }
}

Complex SM_Phantom_U1::f(double tau)
{
  if (tau<=0.) return Complex(0.,0.);
  if (tau>=1.) return Complex(ATOOLS::sqr(::asin(sqrt(1./tau))),0.);

  double eta=sqrt(1.-tau);
  Complex a(log((1.+eta)/(1.-eta)),-M_PI);
  return -.25*a*a;
}

double SM_Phantom_U1::GetFermionContribution(double massratio)
{
  if (massratio<=0.) return 2./3.;
  double tau=ATOOLS::sqr(2.*massratio);
  return tau*(1.+(1.-tau)*real(f(tau)));
}
