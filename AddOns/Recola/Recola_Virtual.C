#include "Recola_Virtual.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace Recola {
 
 Recola_Virtual::Recola_Virtual(const Process_Info& pi,
				const Flavour_Vector& flavs,
				unsigned int recola_id) :
   Virtual_ME2_Base(pi, flavs), m_recola_id(recola_id)
  {
    m_procmap[m_recola_id]=pi;
    m_eventcount=0;
    m_fixedIRscale=true;
    Data_Reader reader(" ",";","#","=");
    m_IRscale=reader.GetValue<int>("IR_SCALE",100);
    m_UVscale=reader.GetValue<int>("UV_SCALE",100);
  }
   
   void Recola_Virtual::Calc(const Vec4D_Vector& momenta) {
    if (Recola_Interface::checkProcGeneration() == false){
     Data_Reader reader(" ",";","#","=");
     int ewscheme = reader.GetValue<int>("EW_SCHEME",1);
     if(ewscheme == 3){
      use_gfermi_scheme_and_set_alpha_rcl(AlphaQED());
     }
     else if (ewscheme == 2){
       //use_alphaz_scheme_and_set_alpha_rcl(AlphaQED());
       use_alphaz_scheme_rcl(AlphaQED());
     }
     else if(ewscheme == 1){
       //use_alpha0_scheme_and_set_alpha_rcl(AlphaQED());
       use_alpha0_scheme_rcl(AlphaQED());
     }
     else{
      cout << "The EW scheme "<<ewscheme<<" is not available with the SHERPA+Recola interface. Valid options are:" << endl;
      cout <<"   1. alpha_QED(0)"<<endl;
      cout <<"   2. alpha_QED(M_Z)"<<endl;
      cout <<"   3. GFermi"<<endl;
      exit(1);
     }
     
     int nlight=0;
     set_mu_ir_rcl(m_IRscale);
     set_mu_uv_rcl(m_UVscale);
     double default_alphas(Recola_Interface::GetDefaultAlphaS());
     double default_scale(Recola_Interface::GetDefaultScale());
     double default_flavscheme(Recola_Interface::GetDefaultFlav());
     double fixed(Recola_Interface::GetFixed());
     double alpha_mat;
     if (default_flavscheme==0)
      default_flavscheme=fixed;
     if (default_flavscheme==16) default_flavscheme=-1;
     if (fixed>0 && fixed<10){
      nlight=fixed;
     }
     else{
      if (default_flavscheme>10){
       nlight=Recola_Interface::PDFnf(m_mur2,default_flavscheme-10);
      }
      if (default_flavscheme==-1)
       nlight=-1;
      if (default_flavscheme==-2 || default_flavscheme==0){
       if (Flavour(kf_c).Mass()!=0)
	nlight=3;
       else if (Flavour(kf_b).Mass()!=0)
	nlight=4;
       else if (Flavour(kf_t).Mass()!=0)
	nlight=5;
       else {
	msg_Out()<<"WARNING: 6 light flavours detected.\n";
	nlight=6;
       }
      }
     }
     if (nlight==0){
      msg_Error()<<METHOD<<"(): Cannot determine number of flavours\n";
     }
     if (nlight>6){
      msg_Error()<<METHOD<<"(): Too many light flavours: "<<nlight<<"\n   Max is 6\n";
     }
     
     if (default_scale==-1. || default_alphas==-1. || Recola_Interface::GetPDFDefault()==0){
      set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),nlight);
      msg_Debugging() << "use AlphaQCD" << AlphaQCD() << "   sqrt(m_mur2) " << sqrt(m_mur2) << "\n";
     }
     else{	
      set_alphas_rcl(default_alphas,sqrt(default_scale),nlight);
      msg_Debugging() << "use default_alphas\n";
     }
     msg_Out() << "processes in Recola are being generated..." << endl;
     generate_processes_rcl();
     Recola_Interface::setProcGenerationTrue();
     msg_Out() << "process generation in Recola completed..." << endl;
     //get_alpha_rcl(alpha_mat);
     //cout << "alpha is: " << 1./alpha_mat << endl;	
     Recola_Interface::SetDefaultFlav(nlight);
    }
    double alpha(0);
    get_alpha_rcl(alpha);
    
    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
     timing = new MyTiming();
     timing->Start();
    }
    set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),Recola_Interface::GetDefaultFlav());
    Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res);
    
    if (msg_LevelIsDebugging()) {
     timing->Stop();
     PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
     <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    }
    
    
    Data_Reader reader(" ",";","#","=");
    m_eventcount+=1;
    
    // factor which by Sherpa convention has to be divided out at this stage
    if (m_born==0) m_born=1.;
    double factor=m_born*AlphaQCD()/2.0/M_PI;
    m_res.Finite()/=factor;
    m_res.IR()/=factor;
    m_res.IR2()/=factor;
    int npoint=reader.GetValue<int>("JUST_ONE_POINT",0);
    if (npoint){
     std::cout<<"QED: "<<AlphaQED()<<std::endl;
     std::cout<<"and scale is:  "<<sqrt(m_mur2)<<std::endl;
     std::cout<<std::setprecision(15)<<"Momenta are:  "<<momenta<<std::endl;
     std::cout<<"flavour:  "<<m_flavs<<std::endl;
     std::cout<< std::setprecision(15)<<"    and finite is:  "<<m_res.Finite()<<std::endl;
     std::cout<< std::setprecision(15)<<"    and finite*factor is:  "<<m_res.Finite()*factor<<std::endl;
     std::cout<< std::setprecision(15)<<"   factor is:  "<<factor<<std::endl;
     std::cout<< std::setprecision(15)<<"   coupling:  "<<AlphaQCD()<<std::endl;
     std::cout<<"And the Born is:  "<<m_born<<std::endl;
     std::cout<<"mu:  "<<sqrt(m_mur2)<<std::endl;
     std::cout<<"AlphaQED():  "<<AlphaQED()<<std::endl;
     std::cout<<"AlphaQCD():  "<<AlphaQCD()<<std::endl;
     std::cout<<"---------------------------"<<std::endl;
     if (m_eventcount==npoint)
      exit(1);
    }
    
   }

}

using namespace Recola;

DECLARE_VIRTUALME2_GETTER(Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Recola_Virtual>::
operator()(const Process_Info &pi) const
{
 DEBUG_FUNC(pi);
 if (pi.m_loopgenerator!="Recola") return NULL;
 
 if (pi.m_fi.m_nloqcdtype!=nlo_type::loop) return NULL;
 
 int procIndex=Recola_Interface::RegisterProcess(pi, 11);
 
 if (procIndex>0) {
  Flavour_Vector flavs = pi.ExtractFlavours();
  return new Recola_Virtual(pi, flavs, procIndex);
 }
 else {
  return NULL;
 }
}
