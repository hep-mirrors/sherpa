#include "MODEL/Main/THDM.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Standard_Model.H"
#include "MODEL/Main/Spectrum_Generator_Base.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(THDM_Getter,"THDM",Model_Base,Model_Arguments);

Model_Base *THDM_Getter::operator()(const Model_Arguments &args) const
{
  return new THDM(args.m_path,args.m_file,args.m_elementary);
}

void THDM_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Two Higgs Doublet Model\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- properties of the h0,H0,A0 and H+ particle\n"
     <<std::setw(width+10)<<" "<<"- MASS[i] and WIDTH[i] for i=25,35,36,37\n"
     <<std::setw(width+7)<<" "<<"- TAN(BETA) (ratio of v_2 and v_1)\n"
     <<std::setw(width+7)<<" "<<"- ALPHA (the Higgs mixing angle)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


THDM::THDM(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool THDM::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the THDM from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("THDM");

  p_sm->ModelInit(isr);
  p_numbers   = p_sm->ExtractScalarNumbers();
  p_constants = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions = p_sm->ExtractScalarFunctions();
  p_matrices  = p_sm->ExtractComplexMatrices();

  delete p_sm;

  FillSpectrum(isr);

  return true;
}

THDM::~THDM()
{
}

void THDM::ParticleInit() {
  //add Higgses
  //kf_code,mass,width,charge,icharge,strong,spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[35] = new Particle_Info(35,1000.,10.0,0,0,0,0,-1,1,1,1,"H0","H_0");
  s_kftable[36] = new Particle_Info(36,1000.,10.0,0,0,0,0,-1,1,1,1,"A0","A_0");
  s_kftable[37] = new Particle_Info(37,1000.,10.0,3,0,0,0,0,1,1,1,"H+","H^+");

  ReadParticleData();
}

void THDM::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
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


  Flavour h0(kf_h0); 
  h0.SetMass(ScalarConstant("Mh0"));
  h0.SetHadMass(ScalarConstant("Mh0"));
  Flavour H0(kf_H0); 
  H0.SetMass(ScalarConstant("MH0"));
  H0.SetHadMass(ScalarConstant("MH0"));
  Flavour A0(kf_A0); 
  A0.SetMass(ScalarConstant("MA0"));
  A0.SetHadMass(ScalarConstant("MA0"));
  Flavour Hplus(kf_Hplus); 
  Hplus.SetMass(ScalarConstant("MHplus"));
  Hplus.SetHadMass(ScalarConstant("MHplus"));

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

