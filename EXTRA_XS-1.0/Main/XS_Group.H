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
    XS_Selector            * p_xsselector;
    bool                     m_atoms;
    std::vector<XS_Base *>   m_xsecs;
    bool                     m_channels;
  public :
    XS_Group(int,int,ATOOLS::Flavour *,
	     PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,ATOOLS::Selector_Data *,
	     int,int,double);
    XS_Group(int,int,ATOOLS::Flavour *);
    XS_Group(int,int,std::string);
    XS_Group();
    ~XS_Group();

    void      PushBack(XS_Base *_xsec) { m_xsecs.push_back(_xsec); }
    void      Add(XS_Base *);
    bool      Delete(XS_Base *);
    void      Clear();
    int       Size() { return m_xsecs.size(); }
    XS_Base * operator[] (int idx) { 
      if ((idx<(int)m_xsecs.size()) && (idx>=0)) return m_xsecs[idx];
      return 0;
    } 
    XS_Selector     * XSSelector() { return p_xsselector; }
    virtual XS_Base * Selected();
    virtual bool      CalculateTotalXSec();
    virtual bool      OneEvent();
    virtual void      SetAtoms(bool _atoms) { m_atoms = _atoms; }
    virtual void      CreateFSRChannels() {}

    void      SetTotalXS();
    void      SetMax(double,int);
    void      AddPoint(const double);
    double    Differential2();
    double    Differential(double,double,double);
    double    operator()(double,double,double)    { return 0.; }
    bool      SetColours(double, double, double)  { return 1; }
    void      UpdateCuts(double,double)           {}
    void      DeSelect();
    void      SelectOne();
    void      SetISR(PDF::ISR_Handler *);
  };
}

#endif


