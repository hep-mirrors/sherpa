#ifndef AMEGIC_Main_Amegic_Base_H
#define AMEGIC_Main_Amegic_Base_H

#include "AMEGIC++/Main/Process_Tags.H"
#include "PHASIC++/Process/Process_Info.H"


namespace AMEGIC {

  class Amegic_Base {
  protected:
    Process_Tags *p_pinfo;
  public:
    Amegic_Base();
    virtual ~Amegic_Base();
    virtual bool PerformTests() = 0;
    virtual bool NewLibs() = 0;
    virtual void Minimize() = 0;

    Pol_Info ExtractPolInfo(const PHASIC::Subprocess_Info &spi);

    void TranslateDecay(Process_Tags &info,const PHASIC::Subprocess_Info &spi);
    Process_Tags *Translate(const PHASIC::Process_Info &pi);
  };
}

#endif
