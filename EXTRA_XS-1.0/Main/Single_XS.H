#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"
#include "Standard_Selector.H"

namespace AMEGIC { class Single_Process; }

namespace EXTRAXS {

  class Single_XS : public XS_Base {
  private:

    AMEGIC::Single_Process * broker;

    double               thres;
    virtual void GenerateName();

  public:

    Single_XS(int,int,APHYTOOLS::Flavour *);
    ~Single_XS(){};

    double            SetThres(double _thres) { return thres = _thres; }
    double            Thres()                 { return thres;          } 

    virtual double    operator() (AMATOOLS::Vec4D *);
    virtual double    operator() (double,double,double);

    virtual XS_Base * Selected() { return this; }

    void MakeBroker(ISR::ISR_Handler *, BEAM::Beam_Spectra_Handler *,
		    APHYTOOLS::Selector_Data *,AMEGIC::Process_Group *);

  };
}

#endif



