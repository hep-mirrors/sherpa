#include "COMIX/Main/Process_Base.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Channels/Sample_Multi_Channel.H"
#include "PHASIC++/Channels/PreSample_Multi_Channel.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Channels/VHAAG.H"
#include "COMIX/Phasespace/PS_Channel.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

std::string COMIX::ComixLogo()
{
  if (!msg->Modifiable()) return "Comix";
  return "\033[31mC\033[32mo\033[34mm\033[0mi\033[33mx\033[0m";
}

COMIX::Process_Base::Process_Base(PHASIC::Process_Base *const proc,Model *const model):
  p_proc(proc), p_model(model), p_psgen(NULL), m_psmc(false),
  m_cls(-1), m_hls(-1), p_pmap(NULL), p_umprocs(NULL) {}

COMIX::Process_Base::~Process_Base() 
{
}

bool COMIX::Process_Base::Initialize(std::map<std::string,std::string> *const pmap,
			   std::vector<Single_Process*> *const procs)
{
  p_pmap=pmap;
  p_umprocs=procs;
  p_proc->Integrator()->SetColorScheme(cls::sample);
  return true;
}

void COMIX::Process_Base::InitModel(MODEL::Model_Base *const model,
			  const std::string &file)
{
  p_model = Model_Getter::GetObject(model->Name(),"");
  if (p_model==NULL) {
    msg_Info()<<METHOD<<"(): No model '"+model->Name()+"' in Comix."<<std::endl;
    return;
  }
  p_model->Initialize(model,file);
}

void COMIX::Process_Base::UpdateIntegrator(Phase_Space_Handler *const psh)
{
  if (!(m_imode&1)) return;
  Multi_Channel *mc(psh->FSRIntegrator());
  if (p_proc->NOut()<=2 || (!m_psmc && mc!=NULL)) return;
  if (m_cls<0) {
    Data_Reader read(" ",";","!","=");
    if (!read.ReadFromFile(m_cls,"CS_OTF")) m_cls=0;
    else msg_Info()<<METHOD<<"(): Set on the flight sampling "<<m_cls<<".\n";
  }
  if (m_hls<0) {
    Data_Reader read(" ",";","!","=");
    if (!read.ReadFromFile(m_hls,"HS_PSMC")) m_hls=2;
    else msg_Info()<<METHOD<<"(): Set helicity sampling mode "<<m_hls<<".\n";
  }
  p_proc->Integrator()->ColorIntegrator()->SetOn(true);
  SP(Helicity_Integrator) helint(p_proc->Integrator()->HelicityIntegrator());
  if (helint!=NULL) {
    if (m_hls&2) {
      p_proc->Integrator()->SetHelicityScheme(hls::sample);
      helint->SetOn(true);
    }
    else {
      p_proc->Integrator()->SetHelicityScheme(hls::sum);
      helint->SetOn(false);
    }
  }
  PreSample_Multi_Channel *psmc(dynamic_cast<PreSample_Multi_Channel*>(mc));
  if (psmc!=NULL) { 
    std::string pid(p_proc->Integrator()->ResultPath()+"/MC_"+p_proc->Name());
    ATOOLS::MakeDir(pid.c_str(),0,false); 
    psmc->WriteOut(pid+"/MC_FSR");
  }
  Channel_Vector chs;
  if (psmc!=NULL) chs=psmc->ExtractChannels();
  Sample_Multi_Channel *smc
    (new Sample_Multi_Channel(psh,&*p_proc->Integrator()->ColorIntegrator()));
  smc->Initialize(chs);
  smc->Reset();
  if (mc!=NULL) {
    smc->SetValidN(mc->ValidN());
    smc->SetN(mc->N());
    mc->DropAllChannels();
    delete mc;
  }
  psh->SetFSRIntegrator(smc);
  p_proc->Integrator()->Reset();
}

bool COMIX::Process_Base::FillIntegrator(Phase_Space_Handler *const psh)
{
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(m_imode,"CDXS_IMODE")) m_imode=2;
  else msg_Info()<<METHOD<<"(): Set integration mode "<<m_imode<<".\n";
  bool pureqcd(true);
  for (size_t i(0);i<p_proc->Flavours().size();++i)
    if (!p_proc->Flavours()[i].Strong()) pureqcd=false;
  if (!pureqcd) m_imode&=~1;
  Multi_Channel *mc(psh->FSRIntegrator());
  if (!(m_imode&1)) {
    mc->DropAllChannels();
    PS_Channel *ch(new PS_Channel(p_proc->NIn(),p_proc->NOut(),
				  (Flavour*)&p_proc->Flavours().front(),this));
    InitPSGenerator(0);
    mc->Add(ch);
    return false;
  }
  SP(Color_Integrator) colint(p_proc->Integrator()->ColorIntegrator());
  if (p_proc->NOut()>2) {
    if (!m_psmc) {
      delete mc;
      psh->SetFSRIntegrator(NULL);
      UpdateIntegrator(psh);
      return false;
    }
    Data_Reader read(" ",";","!","=");
    if (!read.ReadFromFile(m_hls,"HS_PSMC")) m_hls=2;
    else msg_Info()<<METHOD<<"(): Set helicity sampling mode "<<m_hls<<".\n";
    SP(Helicity_Integrator) helint(p_proc->Integrator()->HelicityIntegrator());
    if (helint!=NULL) {
      if (m_hls&1) {
	p_proc->Integrator()->SetHelicityScheme(hls::sample);
	helint->SetOn(true);
      }
      else {
	p_proc->Integrator()->SetHelicityScheme(hls::sum);
	helint->SetOn(false);
      }
    }
    VHAAG *first(NULL);
    PreSample_Multi_Channel *psmc
      (new PreSample_Multi_Channel(psh,&*colint));
    for (size_t i(0);i<(p_proc->NIn()+p_proc->NOut())/2;++i) {
      colint->GenerateType(i,true);
      Idx_Vector perm(colint->Orders().front());
      Multi_Channel *smc(new Multi_Channel(ToString(i),i));
      Idx_Vector rperm(perm.size());
      rperm.front()=perm.front();
      for (size_t j(1);j<rperm.size();++j) rperm[j]=perm[perm.size()-j];
      std::vector<size_t> hp(perm.size());
      for (size_t j(0);j<perm.size();++j) hp[j]=perm[j];
      VHAAG *sc(new VHAAG(p_proc->NIn(),p_proc->NOut(),hp,first));
      if (first==NULL) first=sc;
      smc->Add(sc);
      for (size_t j(0);j<perm.size();++j) hp[j]=rperm[j];
      sc = new VHAAG(p_proc->NIn(),p_proc->NOut(),hp,first);
      smc->Add(sc);
      psmc->AddMC(smc);
    }
    mc->DropAllChannels();
    delete mc;
    psh->SetFSRIntegrator(psmc);
    colint->SetOn(false);
    return true;
  }
  if (!read.ReadFromFile(m_cls,"CS_OTF")) m_cls=0;
  else msg_Info()<<METHOD<<"(): Set on the flight sampling "<<m_cls<<".\n";
  if (m_cls==0) colint->Initialize();
  mc->DropAllChannels();
  mc->Add(new S1Channel(p_proc->NIn(),p_proc->NOut(),
			(Flavour*)&p_proc->Flavours().front()));
  mc->Add(new T1Channel(p_proc->NIn(),p_proc->NOut(),
			(Flavour*)&p_proc->Flavours().front()));
  mc->Add(new U1Channel(p_proc->NIn(),p_proc->NOut(),
			(Flavour*)&p_proc->Flavours().front()));
  return false;
}      

