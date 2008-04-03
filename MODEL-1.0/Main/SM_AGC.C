#include "SM_AGC.H"
#include "Standard_Model.H"
#include "Message.H"
#include "Data_Reader.H"

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(SM_AGC_Getter,"SM+AGC",Model_Base,Model_Arguments);

Model_Base *SM_AGC_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_AGC(args.m_path,args.m_file);
}

void SM_AGC_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Standard Model + Anomalous Gauge Couplings";
}

SM_AGC::SM_AGC(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the Standard Model \\w AGC from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+AGC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();
 
  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;
  
  FillSpectrum();
 
}

void SM_AGC::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  
  //Anomalous gauge couplings (hep-ph/0001065)
  p_constants->insert(std::make_pair(std::string("Alpha_4"),
				     p_dataread->GetValue<double>("ALPHA_4_G_4",0.)));
  p_constants->insert(std::make_pair(std::string("Alpha_5"),
				     p_dataread->GetValue<double>("ALPHA_5",0.)));
  //Anomalous gauge couplings (Nucl. Phys. B282 (1987) 253-307)
  p_constants->insert(std::make_pair(std::string("g1_gamma"),
				     p_dataread->GetValue<double>("G1_GAMMA",1.)));
  p_constants->insert(std::make_pair(std::string("kappa_gamma"),
				     p_dataread->GetValue<double>("KAPPA_GAMMA",1.)));
  p_constants->insert(std::make_pair(std::string("lambda_gamma"),
				     p_dataread->GetValue<double>("LAMBDA_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("g4_gamma"),
				     p_dataread->GetValue<double>("G4_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("g5_gamma"),
				     p_dataread->GetValue<double>("G5_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("kappat_gamma"),
				     p_dataread->GetValue<double>("KAPPAT_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("lambdat_gamma"),
				     p_dataread->GetValue<double>("LAMBDAT_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("g1_Z"),
				     p_dataread->GetValue<double>("G1_Z",1.)));
  p_constants->insert(std::make_pair(std::string("kappa_Z"),
				     p_dataread->GetValue<double>("KAPPA_Z",1.)));
  p_constants->insert(std::make_pair(std::string("lambda_Z"),
				     p_dataread->GetValue<double>("LAMBDA_Z",0.)));
  p_constants->insert(std::make_pair(std::string("g4_Z"),
				     p_dataread->GetValue<double>("G4_Z",0.)));
  p_constants->insert(std::make_pair(std::string("g5_Z"),
				     p_dataread->GetValue<double>("G5_Z",0.)));
  p_constants->insert(std::make_pair(std::string("kappat_Z"),
				     p_dataread->GetValue<double>("KAPPAT_Z",0.)));
  p_constants->insert(std::make_pair(std::string("lambdat_Z"),
				     p_dataread->GetValue<double>("LAMBDAT_Z",0.)));
}


