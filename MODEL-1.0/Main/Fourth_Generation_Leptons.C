#include "Fourth_Generation_Leptons.H"
#include "Message.H"
#include "Standard_Model.H"
#include "Spectrum_Generator_Base.H"

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(Fourth_Generation_Leptons_Getter,"SM+4thLF",Model_Base,Model_Arguments);

Model_Base *Fourth_Generation_Leptons_Getter::operator()(const Model_Arguments &args) const
{
  return new Fourth_Generation_Leptons(args.m_path,args.m_file);
}

void Fourth_Generation_Leptons_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Standard Model + 4th Lepton Family"; 
}

Fourth_Generation_Leptons::Fourth_Generation_Leptons(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the Standard Model \\w 4th Lepton Family from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+4thLF");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  FillSpectrum();

  msg_Tracking()<<METHOD
	   <<" : m(tau') = "<<Flavour(kf_tau_prime).PSMass()<<"("<<Flavour(kf_tau_prime).Charge()<<")"
	   <<",  m(nutau') = "<<Flavour(kf_nutau_prime).PSMass()<<"("<<Flavour(kf_nutau_prime).Charge()<<")"
	   <<"."<<std::endl;
}

Fourth_Generation_Leptons::~Fourth_Generation_Leptons()
{
}

void Fourth_Generation_Leptons::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  double massf4  = p_dataread->GetValue<double>("MASS_L",Flavour(kf_tau_prime).PSMass());
  double massnu4 = p_dataread->GetValue<double>("MASS_NU",Flavour(kf_nutau_prime).PSMass());

  p_constants->insert(std::make_pair(std::string("Yukawa_tauprime"),massf4));
  p_constants->insert(std::make_pair(std::string("Yukawa_nutauprime"),massnu4)); 

  p_numbers->insert(std::make_pair(std::string("CHARGE_L4"),    
				   3*p_dataread->GetValue<int>("CHARGE_L",-1)));
  p_numbers->insert(std::make_pair(std::string("CHARGE_NU4"),    
				   3*p_dataread->GetValue<int>("CHARGE_NU",0)));  

  Flavour taup(kf_tau_prime); taup.SetOn(true);
  Flavour nutaup(kf_nutau_prime); nutaup.SetOn(true);

  taup.SetMass(ScalarConstant("Yukawa_tauprime"));
  nutaup.SetMass(ScalarConstant("Yukawa_nutauprime"));

  taup.SetIntCharge(ScalarNumber("CHARGE_L4"));
  nutaup.SetIntCharge(ScalarNumber("CHARGE_NU4"));
}

