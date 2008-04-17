#ifndef LEV_Base_H
#define LEV_Base_H

#include "Flavour.H"
#include "Vector.H"

namespace EXTRAXS {

  class LEV_Base {
  protected:

    ATOOLS::Flavour m_fla, m_flb, m_flc;

    ATOOLS::Vec4D m_k, m_q;

    double m_kt2min, m_kt2max, m_ymin, m_ymax;

    LEV_Base *p_selected;
    
    // constructors
    inline LEV_Base(const ATOOLS::Flavour &fla,
		    const ATOOLS::Flavour &flb,
		    const ATOOLS::Flavour &flc):
      m_fla(fla), m_flb(flb), m_flc(flc), p_selected(this) {}
    inline LEV_Base(const LEV_Base &base):
      m_fla(base.m_fla), m_flb(base.m_flb), 
      m_flc(base.m_flc), p_selected(this) {}
    inline LEV_Base():
      m_fla(kf_none), m_flb(kf_none), 
      m_flc(kf_none), p_selected(this) {}
    
  public:
    
    // destructor
    virtual ~LEV_Base();

    // member functions
    virtual double Value(const ATOOLS::Vec4D &k1,
			 const ATOOLS::Vec4D &q1,
			 const ATOOLS::Vec4D &k2,
			 const ATOOLS::Vec4D &q2) const;
    virtual double MajorValue(const ATOOLS::Vec4D &k1,
			      const ATOOLS::Vec4D &q1,
			      const ATOOLS::Vec4D &k2,
			      const ATOOLS::Vec4D &q2) const;

    virtual double MajorIntegral(const ATOOLS::Flavour &fl);

    virtual bool SelectSplitting(const double &rn);

    virtual void SetKT2Min(const double &kt2min);
    virtual void SetKT2Max(const double &kt2max);

    virtual LEV_Base *Selected() const;

    // inline functions
    inline ATOOLS::Flavour GetA() const { return m_fla; }
    inline ATOOLS::Flavour GetB() const { return m_flb; }
    inline ATOOLS::Flavour GetC() const { return m_flc; }
      
  };// end of class LEV_Base

}// end of namespace EXTRAXS

#endif
