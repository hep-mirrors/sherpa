#include "MODEL/Main/Hidden_Valley.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Running_Alpha_HV.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(HiddenValley_Getter,"SM+HiddenValley",Model_Base,Model_Arguments);

Model_Base *HiddenValley_Getter::operator()(const Model_Arguments &args) const
{
  return new HiddenValley(args.m_path,args.m_file,args.m_elementary);
}

void HiddenValley_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + Hidden Sector\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- HV_GROUP     (SU, SO or SP)\n"
     <<std::setw(width+7)<<" "<<"- HV_GROUP_DIM (integer value of Nc)\n"
     <<std::setw(width+7)<<" "<<"- ALPHA_HV(MZ) (HV coupling at MZ)\n"
     <<std::setw(width+7)<<" "<<"- ORDER_ALPHA_HV (0: beta0, 1: beta1)\n"
     <<std::setw(width+7)<<" "<<"- MASS[9900001]\n"
     <<std::setw(width+7)<<" "<<"- MASS[9900002]\n"
     <<std::setw(width+7)<<" "<<"- MASS[9900023]\n"
     <<std::setw(width+7)<<" "<<"- MASS[9900024]\n"
     <<std::setw(width+7)<<" "<<"- WIDTH[9900023]\n"
     <<std::setw(width+7)<<" "<<"- WIDTH[9900024]\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


HiddenValley::HiddenValley(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool HiddenValley::ModelInit()
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model Hidden Sector "
	      <<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+HiddenValley");

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

HiddenValley::~HiddenValley() 
{ }

void HiddenValley::ParticleInit() {
  //add two dark flavours, messenger Zp & dark gluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[9900001] = new Particle_Info(9900001,5.,.1,-1,0,3,1,0,1,1,1,"Dark_d","Dark_d");
  s_kftable[9900002] = new Particle_Info(9900002,10.,.1,2,0,3,1,0,1,1,1,"Dark_u","Dark_u");
  s_kftable[9900021] = new Particle_Info(9900021,0.,0.,0,0,8,2,-1,1,0,1,"Dark_g","Dark_g");
  s_kftable[9900023] = new Particle_Info(9900023,200.,10.,0,0,0,2,-1,1,1,1,"Zp","Zp");
  s_kftable[9900024] = new Particle_Info(9900024,200.,10.,3,0,0,2,0,1,1,1,"Wp","Wp");
  ReadParticleData();
}

void HiddenValley::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);

  int    order_alpha_HV = p_dataread->GetValue<int>("ORDER_ALPHA_HV",1);
  double alpha_HV       = p_dataread->GetValue<double>("ALPHA_HV(MZ)",0.118);
  double alpha_HV_def   = p_dataread->GetValue<double>("ALPHA_HV(default)",alpha_HV);
  double MZ2            = sqr((*p_constants)[std::string("MZ")]);
  std::string group     = p_dataread->GetValue<std::string>("HV_GROUP","SU");  
  double HV_Nc          = p_dataread->GetValue<double>("HV_GROUP_DIM",3.);  

  double TR=0.,CA=0.,CF=0.;
  
  if (group==std::string("SU")) {
    CF = (HV_Nc*HV_Nc-1.)/(2.*HV_Nc);        
    CA = HV_Nc;           
    TR = 1./2.;
  }
  else if (group==std::string("SO")) {
    CF = 1./2.*(HV_Nc-1);        
    CA = (HV_Nc-2.);        
    TR=1;
    if (HV_Nc==3) {
      CF *= 2.;
      CA *= 2.;
      TR *= 2.;
    }
  }
  else if (group==std::string("SP")) {
    CF = 1./4.*(HV_Nc+1);
    CA = 0.5*(HV_Nc+2.);     
    TR = 1./2.;
  }
  else {
    msg_Error()<<" Gauge Group not supported assume SU(3) instead! "<<std::endl; 
    HV_Nc=3.;
    CF = (HV_Nc*HV_Nc-1.)/(2.*HV_Nc);        
    CA = HV_Nc;           
    TR = 0.5;
  }
  
  as_HV = new Running_Alpha_HV(alpha_HV,MZ2,order_alpha_HV,1.,group,HV_Nc);
  as_HV->SetDefault(alpha_HV_def);

  msg_Info()<<" ======== Initialize HV gauge coupling ========"<<std::endl; 
  msg_Info()<<"   consider : "<<group<<"("<<HV_Nc<<")"<<" where TR="<<TR<<" CA="<<CA<<" CF="<<CF<< std::endl;  

  p_constants->insert(std::make_pair(std::string("HV_NC"),HV_Nc));
  p_constants->insert(std::make_pair(std::string("TR"),TR));
  p_constants->insert(std::make_pair(std::string("CA"),CA));
  p_constants->insert(std::make_pair(std::string("CF"),CF));

  p_constants->insert(std::make_pair(std::string("alpha_HV(MZ)"),alpha_HV));
  p_functions->insert(std::make_pair(std::string("alpha_HV"),as_HV));

  msg_Info()<<"   set alpha_HV(MZ) to "
	    <<ScalarFunction(string("alpha_HV"),MZ2)<<endl;

   msg_Info()<<"   yields alpha_HV @ 5 GeV "
	    <<ScalarFunction(string("alpha_HV"),sqr(5.))<<endl;
  msg_Info()<<"   yields alpha_HV @ 500 GeV "
	    <<ScalarFunction(string("alpha_HV"),sqr(500.))<<endl;
  msg_Info()<<" =============================================="<<std::endl; 
}



