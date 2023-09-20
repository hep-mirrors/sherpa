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


Recola::Recola_Virtual::Recola_Virtual(const Process_Info& pi,
        const Flavour_Vector& flavs,
        unsigned int recola_id) :
  Virtual_ME2_Base(pi, flavs), m_recola_id(recola_id),
  m_modebackup(m_mode), m_ismapped(false)
{
  m_procmap[m_recola_id]=pi;
  Settings& s = Settings::GetMainSettings();
  m_providespoles=false;
  m_fixedIRscale=true;

  m_IRscale=s["RECOLA_IR_SCALE"].Get<double>();
  m_UVscale=s["RECOLA_UV_SCALE"].Get<double>();
  m_modebackup=m_mode=Recola_Interface::s_vmode;
  m_voqcd = pi.m_maxcpl[0];
  m_boqcd = pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0];

  // init associated contribs
  size_t n(0);
  if (pi.m_fi.m_asscontribs&asscontrib::EW) {
    ++n;
    if (pi.m_fi.m_asscontribs&asscontrib::LO1) {
      ++n;
      if (pi.m_fi.m_asscontribs&asscontrib::LO2) {
        ++n;
        if (pi.m_fi.m_asscontribs&asscontrib::LO3) {
          ++n;
        }
      }
    }
  }
  m_asscontribs.resize(n);  

  if (m_asscontribs.size()>0 && m_voqcd!=m_boqcd+1)
    THROW(fatal_error,"Associated contribs only implemented for NLO QCD.");

}
  
void Recola::Recola_Virtual::Calc(const Vec4D_Vector& momenta) {
  m_mode=m_modebackup;
  if (!Recola_Interface::checkProcGeneration()){
      Recola_Interface::GenerateProcesses(AlphaQED(),AlphaQCD(),
                                          m_IRscale,m_UVscale,m_mur2);
  }

  m_res*=0.; m_born=0.;
  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]=0.;

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }

  double aqcd=AlphaQCD(); 
  int flav=Recola_Interface::GetDefaultFlav();
  set_alphas_rcl(aqcd,sqrt(m_mur2),flav);
  Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res, m_asscontribs);

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

  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]/=factor;
  msg_Debugging()<<"V/B="<<m_res.Finite()<<std::endl;
  for (size_t i(0);i<m_asscontribs.size();++i)
  msg_Debugging()<<"ASS/B="<<m_asscontribs[i]<<std::endl;

}


bool Recola::Recola_Virtual::IsMappableTo(const PHASIC::Process_Info& pi){
  return false;
}

using namespace Recola;

DECLARE_VIRTUALME2_GETTER(Recola::Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,Recola::Recola_Virtual>::
operator()(const PHASIC::Process_Info &pi) const
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
