#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Selector.H"
#include "XS_Base.H"
#include "Selector.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <vector>

namespace EXTRAXS   { class XS_Selector;   }
namespace AMEGIC    { class Process_Group; }
namespace ISR       { class ISR_Handler;   }
namespace BEAM      { class Beam_Handler;  }

using namespace APHYTOOLS;

namespace EXTRAXS {
  class XS_Group : public XS_Base {
  protected:

    AMEGIC::Process_Group * broker;
    APHYTOOLS::Selector_Data * seldata;

    std::vector<XS_Base *> xsecs;
    XS_Selector * xsselector;

  public :

    XS_Group(int,int,std::string);
    XS_Group() {};
    ~XS_Group();

    virtual void       Add(XS_Base *,bool);
    virtual int        Size() { return xsecs.size(); }
    virtual XS_Base *  operator[] (int idx) { 
      if ((idx<xsecs.size()) && (idx>=0)) return xsecs[idx];
      return 0;
    } 
    
    void MakeBroker(ISR::ISR_Handler *, BEAM::Beam_Handler *,
		    APHYTOOLS::Selector_Data *,AMEGIC::Process_Group *);
    AMEGIC::Process_Group * Broker() { return broker; }

  };
}

#endif


