#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"
#include "Standard_Selector.H"
#include "Single_Broker.H"

namespace PHASIC  { class Single_Broker; }

namespace EXTRAXS {

  class Single_XS : public XS_Base {
  protected:
    
    PHASIC::Single_Broker * broker;

    double               thres;

    virtual void GenerateName();

  public:

    Single_XS(int,int,APHYTOOLS::Flavour *);
    ~Single_XS();

    AMEGIC::Process_Base * GetBroker()    { return broker; }
    AMEGIC::Process_Base * CreateBroker() { return broker = new PHASIC::Single_Broker(this); }

    virtual double    operator() (double,double,double);
    virtual XS_Base * Selected() { return this; }

  };
}

#endif



