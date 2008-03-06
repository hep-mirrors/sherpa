#include "Fourth_Generation_Leptons.H"
#include "Message.H"
#include "Standard_Model.H"
#include "Spectrum_Generator_Base.H"

using namespace MODEL;
using namespace ATOOLS;


Fourth_Generation_Leptons::Fourth_Generation_Leptons(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the Fourth_Generation_Leptons from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("Fourth_Generation_Leptons");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  ReadInFile();
  FillMasses();

  std::cout<<METHOD
	   <<" : m(tau') = "<<Flavour(kf_tau_prime).PSMass()<<"("<<Flavour(kf_tau_prime).Charge()<<")"
	   <<",  m(nutau') = "<<Flavour(kf_nutau_prime).PSMass()<<"("<<Flavour(kf_nutau_prime).Charge()<<")"
	   <<"."<<std::endl;
}

Fourth_Generation_Leptons::~Fourth_Generation_Leptons()
{
}

void Fourth_Generation_Leptons::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  double massf4  = p_dataread->GetValue<double>("MASS_L",Flavour(kf_tau_prime).PSMass());
  double massnu4 = p_dataread->GetValue<double>("MASS_NU",Flavour(kf_nutau_prime).PSMass());

  p_constants->insert(std::make_pair(std::string("Yukawa_tauprime"),massf4));
  p_constants->insert(std::make_pair(std::string("Yukawa_nutauprime"),massnu4)); 

  p_numbers->insert(std::make_pair(std::string("CHARGE_L4"),    
				   3*p_dataread->GetValue<int>("CHARGE_L",-1)));
  p_numbers->insert(std::make_pair(std::string("CHARGE_NU4"),    
				   3*p_dataread->GetValue<int>("CHARGE_NU",0)));  
}

void Fourth_Generation_Leptons::FillMasses() {
  Flavour taup(kf_tau_prime); taup.SetOn(true);
  Flavour nutaup(kf_nutau_prime); nutaup.SetOn(true);

  taup.SetMass(ScalarConstant("Yukawa_tauprime"));
  nutaup.SetMass(ScalarConstant("Yukawa_nutauprime"));

  taup.SetIntCharge(ScalarNumber("CHARGE_L4"));
  nutaup.SetIntCharge(ScalarNumber("CHARGE_NU4"));
}
