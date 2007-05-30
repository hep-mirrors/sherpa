#ifndef LO_LEV_H
#define LO_LEV_H

#include "LEV_Base.H"

namespace EXTRAXS {

  class G_GG: public LEV_Base {
  public:
    
    // constructor
    G_GG();

    // member functions
    double Value(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;
    double MajorValue(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		      const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;

    double MajorIntegral(const ATOOLS::Flavour &fl);

  };// end of class G_GG

  class Q_GQ: public LEV_Base {
  public:
    
    // constructor
    Q_GQ(const ATOOLS::Flavour &q);

    // member functions
    double Value(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;
    double MajorValue(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		      const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;

    double MajorIntegral(const ATOOLS::Flavour &fl);

  };// end of class Q_GQ

  class Q_QG: public LEV_Base {
  public:
    
    // constructor
    Q_QG(const ATOOLS::Flavour &q);

    // member functions
    double Value(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;
    double MajorValue(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		      const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;

    double MajorIntegral(const ATOOLS::Flavour &fl);

  };// end of class Q_QG

  class G_QQ: public LEV_Base {
  public:
    
    // constructor
    G_QQ(const ATOOLS::Flavour &q);

    // member functions
    double Value(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;
    double MajorValue(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		      const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2) const;

    double MajorIntegral(const ATOOLS::Flavour &fl);

  };// end of class G_QQ

}// end of namespace EXTRAXS

#endif
