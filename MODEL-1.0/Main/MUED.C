#include "MUED.H"
#include "Message.H"
#include "Standard_Model.H"
#include "MUED_Spectrum.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;


MUED::MUED(string _dir,string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the MUED from "<<m_dir<<" / "<<m_file<<endl;
  m_name      = string("MUED");

  Model_Base * SM = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(SM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(SM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(SM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(SM->GetComplexMatrices()));

  delete SM;

  ReadInFile();
}

void MUED::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);

  m_spectrum  = p_dataread->GetValue<int>("GENERATOR_ON",1);
  m_generator = p_dataread->GetValue<string>("GENERATOR",string("Internal"));
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
}


bool MUED::RunSpectrumGenerator() {
  if (m_spectrum) {
    if (m_generator!=string("Internal")) {
      msg_Error()<<"Error in MUED::RunSpectrumGenerator."<<endl
		 <<"   Unknown spectrum generator : "
		 <<m_generator<<" use internal solution."<<endl;
      m_generator=string("Internal");
    }
    if (m_generator==string("Internal")) {
      p_spectrumgenerator = new MUED_Spectrum(p_dataread,this);
      p_spectrumgenerator->Run("MUED_Spectrum_Output.dat");
      p_spectrumgenerator->FillMasses();
      abort();
      return 1;
    }
    return 0;
  }
  return 1;
}
