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
    bool                     m_atoms;
    XS_Selector            * p_xsselector;
    std::vector<XS_Base *>   m_xsecs;
  private :
    virtual void CreateFSRChannels() {}
  public :
    XS_Group(int,int,APHYTOOLS::Flavour *,
	     ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *,APHYTOOLS::Selector_Data *,
	     int,int,double);
    XS_Group(int,int,APHYTOOLS::Flavour *);
    XS_Group(int,int,std::string);
    XS_Group();
    ~XS_Group();

    void      Add(XS_Base *);
    int       Size() { return m_xsecs.size(); }
    XS_Base * operator[] (int idx) { 
      if ((idx<m_xsecs.size()) && (idx>=0)) return m_xsecs[idx];
      return 0;
    } 
    virtual XS_Base * Selected();
    virtual bool      CalculateTotalXSec();
    virtual bool      OneEvent();
    virtual void      SetAtoms(bool _atoms) { m_atoms = _atoms; }

    void      SetTotalXS();
    void      AddPoint(const double);
    double    Differential2();
    double    Differential(double,double,double);
    double    operator()(double,double,double)    {}
    bool      SetColours(double, double, double)  {}
    void      UpdateCuts(double,double)           {}
    void      DeSelect();
    void      SelectOne();
    void      SetISR(ISR::ISR_Handler *);
  };
}

#endif


