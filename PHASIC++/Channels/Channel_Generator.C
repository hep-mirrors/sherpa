#include "PHASIC++/Channels/Channel_Generator.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Channel_Generator
#define PARAMETER_TYPE PHASIC::Channel_Generator_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;

Channel_Generator::~Channel_Generator()
{
}

void Channel_Generator::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  Getter_Function::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"

namespace PHASIC {

  class Default_Channel_Generator: public Channel_Generator {
  public:
    
    Default_Channel_Generator(const Channel_Generator_Key &key):
      Channel_Generator(key) {}

    int GenerateChannels()
    {
      return p_proc->FillIntegrator
	(&*p_proc->Integrator()->PSHandler());
    }

  };// end of class Default_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(Default_Channel_Generator,"Default",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,Default_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new Default_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    Default_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"ME-generator specific channels";
}
