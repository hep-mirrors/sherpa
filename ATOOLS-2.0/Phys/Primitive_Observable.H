#ifndef Primitive_Observable_h
#define Primitive_Observable_h

#include "Primitive_Observable_Base.H"

#include "Jet_Finder.H"

namespace APHYTOOLS {
  class Differential_Jetrate : public Primitive_Observable_Base {  
  public:
    Differential_Jetrate(int _type,double _xmin,double _xmax,int _nbins,
			 APHYTOOLS::Selector_Base * _sel)
    {
      type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
      name  = std::string("dY_cut");
      histo = new AMATOOLS::Histogram(11,xmin,xmax,nbins);
    };

    Differential_Jetrate(Differential_Jetrate * old) {
      type  = old->Type();
      sel   = old->Sel();
      xmin  = old->Xmin();
      xmax  = old->Xmax();
      nbins = old->Nbins();
      name  = old->Name();
      histo = new AMATOOLS::Histogram(old->Histo());
    }

    Differential_Jetrate * GetCopy() {
      return new Differential_Jetrate(this);
    }

    void Evaluate(double); 
    void Evaluate(int,AMATOOLS::vec4d *,APHYTOOLS::Flavour *,double);
    void Evaluate(const APHYTOOLS::Parton_List &,double);
  };

}

#endif
