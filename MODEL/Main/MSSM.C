#include "MODEL/Main/MSSM.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Standard_Model.H"
#include "MODEL/Main/LesHouches_Interface.H"
#include "MODEL/Main/Spectrum_Generator_Base.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(MSSM_Getter,"MSSM",Model_Base,Model_Arguments);

Model_Base *MSSM_Getter::operator()(const Model_Arguments &args) const
{
  return new MSSM(args.m_path,args.m_file,args.m_elementary);
}

void MSSM_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Minimal Supersymmetric Standard Model\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- SLHA_INPUT (name of Les Houches conform spectrum file)\n"
     <<std::setw(width+7)<<" "<<"- LESHOUCHES_WIDTHS (1 read decay widths from SLHA file (default), 0 ignore widths info)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


MSSM::MSSM(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
}

bool MSSM::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary) 
    msg_Info()<<"Initialize the MSSM from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM");

  p_sm->ModelInit(isr);
  p_numbers   = p_sm->ExtractScalarNumbers();
  p_constants = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions = p_sm->ExtractScalarFunctions();
  p_matrices  = p_sm->ExtractComplexMatrices();

  delete p_sm;

  p_constants->insert(std::make_pair(std::string("mT"),Flavour(kf_t).Yuk()));

  FillSpectrum(isr);
  
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
  
  return true;
}

MSSM::~MSSM() {}


void MSSM::ParticleInit() {
  //add SUSY particles
  //kf_code,mass,width,charge,icharge,strong,spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[35] = new Particle_Info(35,1000.,0.0,0,0,0,0,-1,1,1,1,"H0","H_0");
  s_kftable[36] = new Particle_Info(36,1000.,0.0,0,0,0,0,-1,1,1,1,"A0","A_0");
  s_kftable[37] = new Particle_Info(37,1000.,0.0,3,0,0,0,0,1,1,1,"H+","H^+");
  s_kftable[1000024] = new Particle_Info(1000024,1000.,0.0,3,0,0,1,0,1,1,1,"chargino1","\\chi^+_1");
  s_kftable[1000037] = new Particle_Info(1000037,1000.,0.0,3,0,0,1,0,1,1,1,"chargino2","\\chi^+_2");
  s_kftable[1000022] = new Particle_Info(1000022,1000.,0.0,0,0,0,1,1,1,1,1,"neutralino1","\\chi^0_1");
  s_kftable[1000023] = new Particle_Info(1000023,1000.,0.0,0,0,0,1,1,1,1,1,"neutralino2","\\chi^0_2");
  s_kftable[1000025] = new Particle_Info(1000025,1000.,0.0,0,0,0,1,1,1,1,1,"neutralino3","\\chi^0_3");
  s_kftable[1000035] = new Particle_Info(1000035,1000.,0.0,0,0,0,1,1,1,1,1,"neutralino4","\\chi^0_4");
  s_kftable[1000021] = new Particle_Info(1000021,1000.,0.0,0,0,8,1,1,1,1,1,"gluino","\\tilde{g}");
  s_kftable[1000002] = new Particle_Info(1000002,1000.,0.0,2,1,3,0,0,1,1,1,"sup_L","\\tilde{u}_L");
  s_kftable[1000004] = new Particle_Info(1000004,1000.,0.0,2,1,3,0,0,1,1,1,"scharm_L","\\tilde{c}_L");
  s_kftable[1000006] = new Particle_Info(1000006,1000.,0.0,2,1,3,0,0,1,1,1,"stop_1","\\tilde{t}_1");
  s_kftable[2000002] = new Particle_Info(2000002,1000.,0.0,2,1,3,0,0,1,1,1,"sup_R","\\tilde{u}_R");
  s_kftable[2000004] = new Particle_Info(2000004,1000.,0.0,2,1,3,0,0,1,1,1,"scharm_R","\\tilde{c}_R");
  s_kftable[2000006] = new Particle_Info(2000006,1000.,0.0,2,1,3,0,0,1,1,1,"stop_2","\\tilde{t}_2");
  s_kftable[1000001] = new Particle_Info(1000001,1000.,0.0,-1,-1,3,0,0,1,1,1,"sdown_L","\\tilde{d}_L");
  s_kftable[1000003] = new Particle_Info(1000003,1000.,0.0,-1,-1,3,0,0,1,1,1,"sstrange_L","\\tilde{s}_L");
  s_kftable[1000005] = new Particle_Info(1000005,1000.,0.0,-1,-1,3,0,0,1,1,1,"sbottom_1","\\tilde{b}_1");
  s_kftable[2000001] = new Particle_Info(2000001,1000.,0.0,-1,-1,3,0,0,1,1,1,"sdown_R","\\tilde{d}_R");
  s_kftable[2000003] = new Particle_Info(2000003,1000.,0.0,-1,-1,3,0,0,1,1,1,"sstrange_R","\\tilde{s}_R");
  s_kftable[2000005] = new Particle_Info(2000005,1000.,0.0,-1,-1,3,0,0,1,1,1,"sbottom_2","\\tilde{b}_2");
  s_kftable[1000011] = new Particle_Info(1000011,1000.,0.0,-3,-1,0,0,0,1,1,1,"selectron_L","\\tilde{e}_L");
  s_kftable[1000013] = new Particle_Info(1000013,1000.,0.0,-3,-1,0,0,0,1,1,1,"smuon_L","\\tilde{\\mu}_L");
  s_kftable[1000015] = new Particle_Info(1000015,1000.,0.0,-3,-1,0,0,0,1,1,1,"stau_1","\\tilde{\\tau}_1");
  s_kftable[2000011] = new Particle_Info(2000011,1000.,0.0,-3,-1,0,0,0,1,1,1,"selectron_R","\\tilde{e}_R");
  s_kftable[2000013] = new Particle_Info(2000013,1000.,0.0,-3,-1,0,0,0,1,1,1,"smuon_R","\\tilde{\\mu}_R");
  s_kftable[2000015] = new Particle_Info(2000015,1000.,0.0,-3,-1,0,0,0,1,1,1,"stau_2","\\tilde{\\tau}_2");
  s_kftable[1000012] = new Particle_Info(1000012,1000.,0.0,0,1,0,0,0,1,1,1,"snu_1","\\tilde{\\nu}_1");
  s_kftable[1000014] = new Particle_Info(1000014,1000.,0.0,0,1,0,0,0,1,1,1,"snu_2","\\tilde{\\nu}_2");
  s_kftable[1000016] = new Particle_Info(1000016,1000.,0.0,0,1,0,0,0,1,1,1,"snu_3","\\tilde{\\nu}_3");

  ReadParticleData(); 
}

void MSSM::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
  RunSpectrumGenerator(isr);
}


void MSSM::RunSpectrumGenerator(const PDF::ISR_Handler_Map& isr) {
  p_spectrumgenerator = new LesHouches_Interface(p_dataread,this,m_dir);
  p_spectrumgenerator->Run(isr);
}
