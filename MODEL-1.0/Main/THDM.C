#include "THDM.H"
#include "Message.H"
#include "Standard_Model.H"
#include "Spectrum_Generator_Base.H"

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(THDM_Getter,"THDM",Model_Base,Model_Arguments);

Model_Base *THDM_Getter::operator()(const Model_Arguments &args) const
{
  return new THDM(args.m_path,args.m_file);
}

void THDM_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Two Higgs Doublet Model"; 
}


THDM::THDM(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the THDM from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("THDM");

  Standard_Model * sm = new Standard_Model(m_dir,m_file);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  ReadInFile();
  FillMasses();
}

THDM::~THDM()
{
}

void THDM::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);
  p_constants->insert(std::make_pair(std::string("tan(beta)"),    
				     p_dataread->GetValue<double>("TAN(BETA)",0.)));
  p_constants->insert(std::make_pair(std::string("alpha"),    
				     p_dataread->GetValue<double>("ALPHA",0.)));
  p_constants->insert(std::make_pair(std::string("Mh0"),    
				     p_dataread->GetValue<double>("Mh0",Flavour(kf_h0).Mass())));
  p_constants->insert(std::make_pair(std::string("MH0"),    
				     p_dataread->GetValue<double>("MH0",Flavour(kf_H0).Mass())));
  p_constants->insert(std::make_pair(std::string("MA0"),    
				     p_dataread->GetValue<double>("MA0",Flavour(kf_A0).Mass())));
  p_constants->insert(std::make_pair(std::string("MHplus"),    
				     p_dataread->GetValue<double>("MHplus",Flavour(kf_Hplus).Mass())));
}

void THDM::FillMasses() {
  Flavour h0(kf_h0); h0.SetMass(ScalarConstant("Mh0"));
  Flavour H0(kf_H0); H0.SetMass(ScalarConstant("MH0"));
  Flavour A0(kf_A0); A0.SetMass(ScalarConstant("MA0"));
  Flavour Hplus(kf_Hplus); Hplus.SetMass(ScalarConstant("MHplus"));

  double alpha(ScalarConstant("alpha")),tanb(ScalarConstant("tan(beta)"));
 
  double sina = ::sin(alpha);
  double cosa = cos(alpha);
  
  CMatrix ZR  = CMatrix(2);
  
  ZR[0][0]    = Complex(-sina,0.);
  ZR[0][1]    = Complex(cosa,0.);
  ZR[1][0]    = Complex(cosa,0.);
  ZR[1][1]    = Complex(sina,0.);
  
  double cosb = sqrt(1./(1.+sqr(tanb)));
  double sinb = cosb*tanb;  

  CMatrix ZH  = CMatrix(2);
  
  ZH[0][0]    = Complex(sinb,0.);
  ZH[0][1]    = Complex(-cosb,0.);
  ZH[1][0]    = Complex(cosb,0.);
  ZH[1][1]    = Complex(sinb,0.);

  msg_Tracking()<<"   ZH is : "<<std::endl;
  for (unsigned int i=0;i<2;++i) {
    for (unsigned int j=0;j<2;++j) {
      msg_Tracking()<<"    "<<ZH[i][j]<<" ";
    }
    msg_Tracking()<<std::endl;
  }
  msg_Tracking()<<"   ZR is : "<<std::endl;
  for (unsigned int i=0;i<2;++i) {
    for (unsigned int j=0;j<2;++j) {
      msg_Tracking()<<"    "<<ZR[i][j]<<" ";
    }
    msg_Tracking()<<std::endl;
  }
   
  p_matrices->insert(std::make_pair(std::string("Z_R"),ZR));
  p_matrices->insert(std::make_pair(std::string("Z_H"),ZH));
}
