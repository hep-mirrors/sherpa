#ifndef PHASIC__Channels__Channel_Generator_H
#define PHASIC__Channels__Channel_Generator_H

#include "ATOOLS/Org/Getter_Function.H"

namespace PHASIC {

  class Process_Base;
  class Multi_Channel;

  struct Channel_Generator_Key {
    Process_Base * p_proc;
    Multi_Channel *p_mc;
    Channel_Generator_Key(Process_Base *const proc,
			  Multi_Channel *const mc):
      p_proc(proc), p_mc(mc) {}
  };// end of struct Channel_Generator_Key

  class Channel_Generator {
  protected:

    Multi_Channel *p_mc;
    Process_Base  *p_proc;

  public:

    typedef ATOOLS::Getter_Function
    <Channel_Generator,Channel_Generator_Key> Getter_Function;

  public:

    inline Channel_Generator(const Channel_Generator_Key &key):
    p_mc(key.p_mc), p_proc(key.p_proc) {}

    virtual ~Channel_Generator();

    virtual int GenerateChannels() = 0;

    static void ShowSyntax(const int mode);

  };// end of class Channel_Generator

}// end of namespace PHASIC

#endif

