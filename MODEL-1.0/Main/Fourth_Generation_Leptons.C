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

  Model_Base * SM = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(SM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(SM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(SM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(SM->GetComplexMatrices()));

  ReadInFile();
  FillMasses();

  std::cout<<METHOD
	   <<" : m(tau') = "<<Flavour(kf::tau_prime).PSMass()<<"("<<Flavour(kf::tau_prime).Charge()<<")"
	   <<",  m(nutau') = "<<Flavour(kf::nutau_prime).PSMass()<<"("<<Flavour(kf::nutau_prime).Charge()<<")"
	   <<"."<<std::endl;
}

void Fourth_Generation_Leptons::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  double massf4  = p_dataread->GetValue<double>("MASS_L",Flavour(kf::tau_prime).PSMass());
  double massnu4 = p_dataread->GetValue<double>("MASS_NU",Flavour(kf::nutau_prime).PSMass());

  p_constants->insert(std::make_pair(std::string("Yukawa_tauprime"),massf4));
  p_constants->insert(std::make_pair(std::string("Yukawa_nutauprime"),massnu4)); 

  p_numbers->insert(std::make_pair(std::string("CHARGE_L4"),    
				   3*p_dataread->GetValue<int>("CHARGE_L",-1)));
  p_numbers->insert(std::make_pair(std::string("CHARGE_NU4"),    
				   3*p_dataread->GetValue<int>("CHARGE_NU",0)));  
}

void Fourth_Generation_Leptons::FillMasses() {
  Flavour(kf::tau_prime).SetOn(true);
  Flavour(kf::nutau_prime).SetOn(true);

  Flavour(kf::tau_prime).SetMass(ScalarConstant("Yukawa_tauprime"));
  Flavour(kf::nutau_prime).SetMass(ScalarConstant("Yukawa_nutauprime"));

  Flavour(kf::tau_prime).SetIntCharge(ScalarNumber("CHARGE_L4"));
  Flavour(kf::nutau_prime).SetIntCharge(ScalarNumber("CHARGE_NU4"));
}
