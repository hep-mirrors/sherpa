#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"
#include "Standard_Selector.H"

namespace EXTRAXS {

  class Single_XS : public XS_Base {
  protected:
    double m_thres;
  public:
    Single_XS(int,int,APHYTOOLS::Flavour *,
	      PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,APHYTOOLS::Selector_Data *,
	      int,int,double);
    Single_XS(int,int,APHYTOOLS::Flavour *);

    virtual int                        Size();
    virtual XS_Base                  * operator[](int idx);
    virtual void                       SetISR(PDF::ISR_Handler *);

    bool                               CalculateTotalXSec();
    void                               SetTotalXS();
    void                               AddPoint(const double);

    virtual double                     Differential2();
    virtual double                     Differential(double,double,double);
    virtual double                     operator()(double,double,double)   = 0;
    virtual bool                       SetColours(double, double, double) = 0;

    virtual bool                       OneEvent();
    virtual void                       UpdateCuts(double,double)           {}
    virtual void                       DeSelect()                          {}
    virtual void                       SelectOne()                         {}
  };

  inline int Single_XS::Size() { return 1; }
  inline XS_Base * Single_XS::operator[](int idx) {
    AORGTOOLS::msg.Error()<<"Error in XS_Base * Single_XS::operator[]. "<<endl
			  <<"   This is absolutely meaningless. Return this."<<endl;
    return this;
  }
  inline void Single_XS::SetISR(PDF::ISR_Handler * _isr) { p_isr     = _isr; }
}

#endif



