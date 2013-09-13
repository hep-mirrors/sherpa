#include "MODEL/Main/SM_AGC.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(SM_AGC_Getter,"SM+AGC",Model_Base,Model_Arguments);

Model_Base *SM_AGC_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_AGC(args.m_path,args.m_file,args.m_elementary);
}

void SM_AGC_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
   str<<"The Standard Model + Anomalous Gauge Couplings\n"
      <<std::setw(width+4)<<" "<<"{\n"
      <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
      <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
      <<std::setw(width+7)<<" "<<"- ALPHA_4_G_4 and ALPHA_5 (L_4,5 parameter of Ref. hep-ph/0001065)\n"
      <<std::setw(width+7)<<" "<<"- Lagrangian parameters of Hagiwara et. al, Nucl. Phys. B282:253,1987\n"
      <<std::setw(width+10)<<" "<<"W-W-Z/photon:\n"
      <<std::setw(width+10)<<" "<<"- G1_GAMMA,KAPPA_GAMMA,LAMBDA_GAMMA,G4_GAMMA,G5_GAMMA\n"
      <<std::setw(width+10)<<" "<<"- KAPPAT_GAMMA,LAMBDAT_GAMMA,G1_Z,KAPPA_Z,LAMBDA_Z\n"
      <<std::setw(width+10)<<" "<<"- G4_Z,G5_Z,KAPPAT_Z,LAMBDAT_Z\n"
      <<std::setw(width+10)<<" "<<"Z-Z/photon-Z/photon:\n"
      <<std::setw(width+10)<<" "<<"- F4_GAMMA,F5_GAMMA,F4_Z,F5_Z\n"
      <<std::setw(width+10)<<" "<<"- H1_GAMMA,H2_GAMMA,H3_GAMMA,H4_GAMMA\n"
      <<std::setw(width+10)<<" "<<"- H1_Z,H2_Z,H3_Z,H4_Z\n"
      <<std::setw(width+7)<<" "<<"- Parameters for unitarity form factor: UNITARIZATION_SCALE, UNITARIZATION_N\n"
      <<std::setw(width+7)<<" "<<"  (see U. Baur, D. Zeppenfeld Nucl.Phys.B308:127,1988)\n"
      <<std::setw(width+4)<<" "<<"}\n";
}

SM_AGC::SM_AGC(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_AGC::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model \\w AGC from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+AGC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();
 
  p_sm->ModelInit(isr);
  p_numbers   = p_sm->ExtractScalarNumbers();
  p_constants = p_sm->ExtractScalarConstants();
  p_functions = p_sm->ExtractScalarFunctions();
  p_matrices  = p_sm->ExtractComplexMatrices();

  delete p_sm;
  
  FillSpectrum(isr);

  return true;
}

void SM_AGC::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
  
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
  p_constants->insert(std::make_pair(std::string("f4_gamma"),
				     p_dataread->GetValue<double>("F4_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("f5_gamma"),
				     p_dataread->GetValue<double>("F5_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("f4_Z"),
				     p_dataread->GetValue<double>("F4_Z",0.)));
  p_constants->insert(std::make_pair(std::string("f5_Z"),
				     p_dataread->GetValue<double>("F5_Z",0.)));
  p_constants->insert(std::make_pair(std::string("h1_gamma"),
				     p_dataread->GetValue<double>("H1_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("h2_gamma"),
				     p_dataread->GetValue<double>("H2_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("h3_gamma"),
				     p_dataread->GetValue<double>("H3_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("h4_gamma"),
				     p_dataread->GetValue<double>("H4_GAMMA",0.)));
  p_constants->insert(std::make_pair(std::string("h1_Z"),
				     p_dataread->GetValue<double>("H1_Z",0.)));
  p_constants->insert(std::make_pair(std::string("h2_Z"),
				     p_dataread->GetValue<double>("H2_Z",0.)));
  p_constants->insert(std::make_pair(std::string("h3_Z"),
				     p_dataread->GetValue<double>("H3_Z",0.)));
  p_constants->insert(std::make_pair(std::string("h4_Z"),
				     p_dataread->GetValue<double>("H4_Z",0.)));
  double n(p_dataread->GetValue<double>("UNITARIZATION_N",2.));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_N"),n));
  double mu(p_dataread->GetValue<double>("UNITARIZATION_SCALE",1000.));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_SCALE"),mu));
  double n3(p_dataread->GetValue<double>("UNITARIZATION_N3",n));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_N3"),n3));
  double mu3(p_dataread->GetValue<double>("UNITARIZATION_SCALE3",mu));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_SCALE3"),mu3));
  double n4(p_dataread->GetValue<double>("UNITARIZATION_N4",n));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_N4"),n4));
  double mu4(p_dataread->GetValue<double>("UNITARIZATION_SCALE4",mu));
  p_constants->insert(std::make_pair(std::string("UNITARIZATION_SCALE4"),mu4));
}


