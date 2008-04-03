#include "SM_Phantom_U1.H"
#include "Standard_Model.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(SM_Phantom_U1_Getter,"PHANTOM_U1",Model_Base,Model_Arguments);

Model_Base *SM_Phantom_U1_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_Phantom_U1(args.m_path,args.m_file);
}

void SM_Phantom_U1_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Standard Model + phantom Higgs"; 
}


SM_Phantom_U1::SM_Phantom_U1(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the SM_Phantom_U1 from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM_Phantom_U1");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  FillSpectrum();

  if (!SanityChecks()) {
    msg_Error()<<"Potential Error in "<<METHOD<<":"<<endl
	       <<"   Sanity checks not passed."<<endl
	       <<"   Continue and hope for the best."<<endl;
  }
}

SM_Phantom_U1::~SM_Phantom_U1() 
{ }

void SM_Phantom_U1::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
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
  HiggsMix[0][1] = sqrt(1.-sqr(abs(HiggsMix[1][1])));
  HiggsMix[1][0] = -HiggsMix[0][1];

  p_matrices->insert(std::make_pair(std::string("HiggsMix"),HiggsMix));
  
  Flavour flav;
  flav = Flavour(kf_h0);
  flav.SetMass(ScalarConstant(string("M_H1")));
  flav.SetMassOn(true);
  flav = Flavour(kf_H0);
  flav.SetMass(ScalarConstant(string("M_H2")));
  flav.SetMassOn(true);
  flav = Flavour(kf_A0);
  flav.SetMass(0.);
  flav.SetMassOn(true);
  flav = Flavour(kf_Z0_2);
  flav.SetMass(ScalarConstant(string("M_Z'")));
  flav.SetMassOn(true);
}

bool SM_Phantom_U1::SanityChecks() {
  if (ScalarConstant(string("Tan(Beta)"))<1.e-6 ||
      ScalarConstant(string("M_H1"))<0.       || 
      ScalarConstant(string("M_H2"))<0.) return false;
  return true;
}



