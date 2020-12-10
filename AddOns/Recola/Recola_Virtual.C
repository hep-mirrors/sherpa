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
    Virtual_ME2_Base(pi, flavs), m_recola_id(recola_id),
    m_modebackup(m_mode), m_ismapped(false)
  {
    //set_mu_ir_rcl(m_mur);
    //set_mu_uv_rcl(m_mur);
    m_procmap[m_recola_id]=pi;
    //m_providespoles=false;
    //m_fixedIRscale=true;
    Settings& s = Settings::GetMainSettings();
    m_IRscale=s["RECOLA_IR_SCALE"].Get<double>();
    m_UVscale=s["RECOLA_UV_SCALE"].Get<double>();
    m_modebackup=m_mode=Recola_Interface::s_vmode;
  }
  
  void Recola_Virtual::Calc(const Vec4D_Vector& momenta) {
    m_mode=m_modebackup;
    if (!Recola_Interface::checkProcGeneration()){
        std::cout<<"process generation started..."<<std::endl;
        Recola_Interface::GenerateProcesses(AlphaQED(),AlphaQCD(),
                                            m_IRscale,m_UVscale,m_mur2);
    }
    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
      timing = new MyTiming();
      timing->Start();
    }

    double aqcd=AlphaQCD(); 
    int flav=Recola_Interface::GetDefaultFlav();
    set_alphas_rcl(aqcd,sqrt(m_mur2),flav);

    Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res);

    if (msg_LevelIsDebugging()) {
      timing->Stop();
      PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
		 <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    }
 
    
    double coupling(1.);
    if (m_stype&sbt::qcd) coupling=AlphaQCD();
    else if (m_stype&sbt::qed) coupling=AlphaQED();
    else THROW(fatal_error,"Unknown coupling.");
    

    // if Born vanishes, do not divide by it, reset mode for this event
    if(!(m_mode&1) && m_born==0.) {
      m_mode|=1;
      msg_Tracking()<<METHOD<<"(): switch to mode 1, Born vanishes"<<std::endl;
    }
    double factor=((m_mode&1)?1.:m_born)*coupling/2.0/M_PI;
    msg_Debugging()<<"cpl="<<coupling/2.0/M_PI<<std::endl;
    // factor which by Sherpa convention has to be divided out at this stage
    m_res.Finite()/=factor;
    m_res.IR()/=factor;
    m_res.IR2()/=factor;

   }


  bool Recola_Virtual::IsMappableTo(const PHASIC::Process_Info& pi){
    Process_Info looppi(pi);
    Process_Info mappi(Recola_Interface::s_procmap[m_recola_id]);
    if (looppi.m_fi.m_nlotype!=nlo_type::lo) looppi.m_fi.m_nlotype=nlo_type::loop;

    std::string nameloop(PHASIC::Process_Base::GenerateName(looppi.m_ii,looppi.m_fi));
    std::string namemap(PHASIC::Process_Base::GenerateName(mappi.m_ii,mappi.m_fi));

    DEBUG_FUNC(nameloop);
    if (namemap==nameloop) m_ismapped=true;
    else                   m_ismapped=false;
    msg_Debugging()<<(m_ismapped?"yes":"no")<<std::endl;

    return m_ismapped;
  }

}

using namespace Recola;

DECLARE_VIRTUALME2_GETTER(Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Recola_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;

  if (pi.m_fi.m_nlotype!=nlo_type::loop) return NULL;

  int procIndex=Recola_Interface::RegisterProcess(pi, 11);
  
  if (procIndex>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Virtual(pi, flavs, procIndex);
  }
  else {
    return NULL;
  }
}
