#include "Fourth_Generation_Leptons.H"
#include "Message.H"
#include "Standard_Model.H"
#include "Spectrum_Generator_Base.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(Fourth_Generation_Leptons_Getter,"SM+4thLF",Model_Base,Model_Arguments);

Model_Base *Fourth_Generation_Leptons_Getter::operator()(const Model_Arguments &args) const
{
  return new Fourth_Generation_Leptons(args.m_path,args.m_file,args.m_elementary);
}

void Fourth_Generation_Leptons_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + 4th Lepton Family\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- properties of \\tau' and \\nu_\\tau' particle\n"
     <<std::setw(width+10)<<" "<<"- MASS[17], MASS[18], WIDTH[17] and WIDTH[18]\n"
     <<std::setw(width+4)<<" "<<"}\n";
    
}

Fourth_Generation_Leptons::Fourth_Generation_Leptons(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model \\w 4th Lepton Family from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+4thLF");
  
  Standard_Model * sm = new Standard_Model(m_dir,m_file,false);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;

  ParticleInit();
  FillSpectrum();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

Fourth_Generation_Leptons::~Fourth_Generation_Leptons()
{
}

void Fourth_Generation_Leptons::ParticleInit() {
  //add 4th lepton generation
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_tau_prime] = new Particle_Info(kf_tau_prime,1000.,10.0,-3,-1,0,1,0,1,1,1,"tau-'","\\tau^{\\prime-}");
  s_kftable[kf_nutau_prime] = new Particle_Info(kf_nutau_prime,1000.,10.0,0,1,0,1,0,1,1,0,"nu_tau'","\\nu_{\\tau'}");
  
  ReadParticleData();  
}

void Fourth_Generation_Leptons::FillSpectrum() {
  double massf4  = Flavour(kf_tau_prime).PSMass();
  double massnu4 = Flavour(kf_nutau_prime).PSMass();

  p_constants->insert(std::make_pair(std::string("Yukawa_tauprime"),massf4));
  p_constants->insert(std::make_pair(std::string("Yukawa_nutauprime"),massnu4)); 
}

