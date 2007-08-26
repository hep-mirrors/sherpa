#include "MSSM_LQQ.H"
#include "MSSM.H"
#include "Message.H"

using namespace MODEL;
using namespace ATOOLS;


MSSM_LQQ::MSSM_LQQ(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the MSSM_LQQ from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM_LQQ");

  Model_Base * MSSM = new MSSM(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(MSSM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(MSSM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(MSSM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(MSSM->GetComplexMatrices()));
  delete MSSM;

  ReadInFile();
}

void MSSM_LQQ::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  CMatrix LQQ_1(3);

  for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) CKM[i][j] = CKM[j][i] = Complex(0.,0.);
    CKM[i][i] = Complex(1.,0.);
  }
  
  double Cabibbo=0.0,A=.8,rho,eta;
  m_ckmorder     = p_dataread->GetValue<int>("CKMORDER",0);  
  if (m_ckmorder>0) {
    Cabibbo    = p_dataread->GetValue<double>("CABIBBO",0.22);
    CKM[0][0] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[1][1] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[0][1] += Cabibbo * Complex( 1.,0.);
    CKM[1][0] += Cabibbo * Complex(-1.,0.);
  }
  if (m_ckmorder>1) {
    A          = p_dataread->GetValue<double>("A",0.8);
    CKM[1][2] += A*sqr(Cabibbo)  * Complex( 1.,0.);
    CKM[2][1] += A*sqr(Cabibbo)  * Complex(-1.,0.);
  }
  if (m_ckmorder>2) {
    eta        = p_dataread->GetValue<double>("ETA",0.5);
    rho        = p_dataread->GetValue<double>("RHO",0.5);
    CKM[0][2] += A*pow(Cabibbo,3) * Complex(rho,-eta);
    CKM[2][0] += A*pow(Cabibbo,3) * Complex(1.-rho,-eta);
  }
  p_matrices->insert(std::make_pair(std::string("CKM"),CKM));
}


bool MSSM_LQQ::RunSpectrumGenerator() {
  if (m_spectrum) {
    m_generator = p_dataread->GetValue<std::string>("SUSY_GENERATOR",std::string("LesHouches"));
#ifdef USING__ISAJET
    if (m_generator==std::string("Isajet")) {
      p_spectrumgenerator = new ISAJET::Isajet_Fortran_Interface(p_dataread,this);
      p_spectrumgenerator->Run(std::string(m_scenario));
      p_spectrumgenerator->FillMasses();
      //p_spectrumgenerator->FillDecays();
      return 1;
    }
#endif
    if (m_generator==std::string("LesHouches")) {
      p_spectrumgenerator = new LesHouches_Interface(p_dataread,this,m_dir);
      p_spectrumgenerator->Run(std::string(m_scenario));
      p_spectrumgenerator->FillMasses();
      //p_spectrumgenerator->FillDecays();
      return 1;
    }
    
    msg_Error()<<"Error in MSSM_LQQ::RunSpectrumGenerator."<<std::endl
	       <<"   Unknown spectrum generator : "<<m_generator<<" use internal solution."<<std::endl;
    return 0;
  }
  return 1;
}
