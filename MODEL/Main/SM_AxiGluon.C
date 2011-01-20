#include "MODEL/Main/SM_AxiGluon.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(SM_AxiGluon_Getter,"SM+AxiGluon",Model_Base,Model_Arguments);

Model_Base *SM_AxiGluon_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_AxiGluon(args.m_path,args.m_file,args.m_elementary);
}

void SM_AxiGluon_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + axigluon\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- MASS_AXIGLUON \n"
     <<std::setw(width+7)<<" "<<"- ...\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


SM_AxiGluon::SM_AxiGluon(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_AxiGluon::ModelInit()
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model plus U(1) phantom Higgs from "
	      <<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+AxiGluon");

  p_sm->ModelInit();
  p_numbers          = p_sm->ExtractScalarNumbers();
  p_constants        = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions        = p_sm->ExtractScalarFunctions();
  p_matrices         = p_sm->ExtractComplexMatrices();

  delete p_sm;

  FillSpectrum();

  return true;
}

SM_AxiGluon::~SM_AxiGluon() 
{ }

void SM_AxiGluon::ParticleInit() {
  //add axigluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[61] = new Particle_Info(61,1000.,10.0,0,0,8,2,-1,1,0,1,"axigluon","\\chi_g");

  ReadParticleData();
}

void SM_AxiGluon::FillSpectrum() {
  p_dataread->RereadInFile();
  p_constants->insert(make_pair(string("MASS_AXI"),    
				p_dataread->GetValue<double>("MASS_AXIGLUON",1000.)));
  Flavour flav = Flavour(61);
  flav.SetMass(ScalarConstant(string("MASS_AXI")));
  flav.SetHadMass(ScalarConstant(string("MASS_AXI")));
  flav.SetMassOn(true);
  flav.SetStable(false);
  flav.SetWidth(-1.);
}
