#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Selector.H"
#include "XS_Base.H"


#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <vector>

namespace EXTRAXS { class XS_Selector; }

namespace EXTRAXS {
  class XS_Group : public XS_Base {
  protected:
    std::vector<XS_Base *> xs;
    XS_Selector * xsselector;
  public :
    /*------------------------------------------------------------------------------
      
      Constructors
      
      ------------------------------------------------------------------------------*/
    XS_Group(int,int,std::string);
    XS_Group() {};
    ~XS_Group();
    /*------------------------------------------------------------------------------
      
      Management of XS_Groups
      
      ------------------------------------------------------------------------------*/
    virtual void Add(XS_Base *);
    virtual int        Size() { return xs.size(); }
    virtual XS_Base *  operator[] (int idx) { 
      if ((idx<xs.size()) && (idx>=0)) return xs[idx];
      return 0;
    } 


    void SetTotalXS();
    void SelectOne();
    void DeSelect();
    XS_Base * Selected();

    virtual void   AddPoint(double);
    virtual double Differential(AMATOOLS::vec4d *);
    virtual double Differential2();
    virtual double Differential(double,double,double);
    virtual double DSigma(double,double,double);
    virtual double DSigma2();


    virtual void SetISR(ISR::ISR_Handler * _isr);
  };
}

#endif


