#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"

using namespace PHASIC;
using namespace ATOOLS;

#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)
typedef PHASIC::Single_Channel *(*Lib_Getter_Function)
  (int nin,int nout,ATOOLS::Flavour* fl,
   ATOOLS::Integration_Info * const info,PHASIC::Phase_Space_Handler *psh);

namespace PHASIC {
  class fsrchannels2_Channel_Generator: public Channel_Generator {
  public:
    fsrchannels2_Channel_Generator(const Channel_Generator_Key &key):
      Channel_Generator(key) {}
    Single_Channel *LoadChannel(int nin,int nout,Flavour* fl,const std::string &pID,Phase_Space_Handler *psh)
    {
      size_t pos(pID.find("/"));
      s_loader->AddPath(rpa->gen.Variable("SHERPA_LIB_PATH"));
      Lib_Getter_Function gf = (Lib_Getter_Function)
        PT(s_loader->GetLibraryFunction("SherpaHiggs","Getter_"+pID));
      if (gf==NULL) return NULL;
      return gf(nin,nout,fl,psh->GetInfo(),psh);
    }
    int GenerateChannels()
    {
      Data_Reader read(" ",";","#","=");
      int os=read.GetValue<int>("HIGGS_ON_SHELL",0);
      int nin=p_proc->NIn(), nout=p_proc->NOut();
      Flavour *fl=(Flavour*)&p_proc->Flavours().front();
      Phase_Space_Handler *psh=&*p_proc->Integrator()->PSHandler();
      p_mc->Add(LoadChannel(nin,nout,fl,"C2_2",psh));
      if (!os) {
      p_mc->Add(LoadChannel(nin,nout,fl,"C2_1",psh));
      p_mc->Add(LoadChannel(nin,nout,fl,"C2_0",psh));
      }
      return 0;
    }
  };
}

DECLARE_GETTER(fsrchannels2_Channel_Generator,"PS2",Channel_Generator,Channel_Generator_Key);
Channel_Generator *ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,fsrchannels2_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const { return new fsrchannels2_Channel_Generator(args); }
void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,fsrchannels2_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"PS2"; }


