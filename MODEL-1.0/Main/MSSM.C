#include "MSSM.H"
#include "Message.H"
#include "Standard_Model.H"
#include "LesHouches_Interface.H"
#include "Spectrum_Generator_Base.H"

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(MSSM_Getter,"MSSM",Model_Base,Model_Arguments);

Model_Base *MSSM_Getter::operator()(const Model_Arguments &args) const
{
  return new MSSM(args.m_path,args.m_file);
}

void MSSM_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MSSM"; 
}


MSSM::MSSM(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the MSSM from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  p_constants->insert(std::make_pair(std::string("mT"),    
				     ScalarConstant("Yukawa_t")));

  FillSpectrum();
}

MSSM::~MSSM() {}

void MSSM::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  RunSpectrumGenerator();
}


void MSSM::RunSpectrumGenerator() {
  p_spectrumgenerator = new LesHouches_Interface(p_dataread,this,m_dir);
  p_spectrumgenerator->Run();
}
