#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Selector.H"
#include "XS_Base.H"
#include "Broker_Group.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <vector>

namespace EXTRAXS { class XS_Selector;  }
namespace PHASIC  { class Broker_Group; }

namespace EXTRAXS {
  class XS_Group : public XS_Base {
  protected:

    std::vector<XS_Base *> xsecs;
    XS_Selector * xsselector;

    PHASIC::Broker_Group * broker;

  public :

    XS_Group(int,int,std::string);
    XS_Group() {};
    ~XS_Group();

    AMEGIC::Process_Base * GetBroker() { return broker; }
    AMEGIC::Process_Base * CreateBroker();

    virtual void       Add(XS_Base *);
    virtual int        Size() { return xsecs.size(); }
    virtual XS_Base *  operator[] (int idx) { 
      if ((idx<xsecs.size()) && (idx>=0)) return xsecs[idx];
      return 0;
    } 
  };
}

#endif


