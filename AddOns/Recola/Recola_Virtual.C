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
  }
  
  void Recola_Virtual::Calc(const Vec4D_Vector& momenta) {
    if (Recola_Interface::checkProcGeneration() == false){
      generate_processes_rcl();
      Recola_Interface::setProcGenerationTrue();
    }
    double alpha(0);
    get_alpha_rcl(alpha);

    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
      timing = new MyTiming();
      timing->Start();
    }

    double aqcd=AlphaQCD(); double mur=sqrt(m_mur2); int flav=Recola_Interface::GetDefaultFlav();
    set_alphas_rcl(aqcd,mur,flav);

    Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res);

    if (msg_LevelIsDebugging()) {
      timing->Stop();
      PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
		 <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    }


    Data_Reader reader(" ",";","#","=");
    m_eventcount+=1;
    
    double coupling(1.);
    if (m_stype&sbt::qcd) {
      coupling=AlphaQCD();
    }
    else if (m_stype&sbt::qed) {
      coupling=alpha;
    }
    else THROW(fatal_error,"Unknown coupling.");
    // factor which by Sherpa convention has to be divided out at this stage
    if (m_born==0) m_born=1.;
    double factor=m_born*coupling/2.0/M_PI;
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
      std::cout<< std::setprecision(15)<<"   coupling:  "<<coupling<<std::endl;
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
