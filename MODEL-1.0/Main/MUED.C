#include "MUED.H"
#include "Message.H"
#include "Standard_Model.H"
#include "MUED_Spectrum.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;


DECLARE_GETTER(MUED_Getter,"MUED",Model_Base,Model_Arguments);

Model_Base *MUED_Getter::operator()(const Model_Arguments &args) const
{
  return new MUED(args.m_path,args.m_file);
}

void MUED_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MUED"; 
}

MUED::MUED(string _dir,string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the MUED from "<<m_dir<<" / "<<m_file<<endl;
  m_name      = string("MUED");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  FillSpectrum();
}

MUED::~MUED()
{
}

void MUED::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);

  p_numbers->insert(make_pair(string("N_Generations"),
			      p_dataread->GetValue<int>("GENERATIONS",1)));
  p_constants->insert(make_pair(string("1/Radius"),
				p_dataread->GetValue<double>("1/RADIUS",1000.)));
  p_constants->insert(make_pair(string("Lambda"),
				p_dataread->GetValue<double>("CUTOFF",20000.)));
  p_constants->insert(make_pair(string("M2bar_H"),
				p_dataread->GetValue<double>("M2_H(BAR)",0.)));

  cout<<METHOD<<": Test this : "
      <<ScalarNumber("N_Generations")<<" / "
      <<ScalarConstant("1/Radius")<<" / "
      <<ScalarConstant("Lambda")<<"."<<endl;

  RunSpectrumGenerator();
}

void MUED::RunSpectrumGenerator() {
  p_spectrumgenerator = new MUED_Spectrum(p_dataread,this);
  p_spectrumgenerator->Run();
}
