#ifndef Primitive_Observable_h
#define Primitive_Observable_h

#include "Primitive_Observable_Base.H"
#include "Jet_Finder.H"

namespace ATOOLS {
  class Differential_Jetrate : public Primitive_Observable_Base {  
  public:
    Differential_Jetrate(int _type,double _xmin,double _xmax,int _nbins,
			 Selector_Base * _sel)
    {
      m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; p_sel = _sel;
      m_name  = std::string("dY_cut");
      p_histo = new Histogram(11,m_xmin,m_xmax,m_nbins);
    };

    Differential_Jetrate(Differential_Jetrate * old) {
      m_type  = old->Type();
      p_sel   = old->Sel();
      m_xmin  = old->Xmin();
      m_xmax  = old->Xmax();
      m_nbins = old->Nbins();
      m_name  = old->Name();
      p_histo = new Histogram(old->Histo());
    }

    Differential_Jetrate * Copy() const {
      return new Differential_Jetrate(m_type,m_xmin,m_xmax,m_nbins,p_sel);
    }

    void Evaluate(int,const Vec4D *,const Flavour *,double, int);
    void Evaluate(const Particle_List &,double, int);
  };

}

#endif
