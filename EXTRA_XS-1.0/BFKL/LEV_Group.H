#ifndef LEV_Group_H
#define LEV_Group_H

#include "LEV_Base.H"

namespace EXTRAXS {

  class LEV_Group: public LEV_Base {
  public:

    typedef std::vector<LEV_Base*> LEV_Vector;
    typedef std::vector<double>    Double_Vector;

  protected:

    LEV_Vector    m_levs;
    Double_Vector m_integrals;

  public:
    
    // constructor
    LEV_Group();

    // destructor
    ~LEV_Group();

    // member functions
    double Value(const double &y,const double &kt) const;
    double MajorValue(const double &y,const double &kt) const;

    double MajorIntegral();

    void SetYMin(const double &ymin);
    void SetYMax(const double &ymax);
    
    void SetKT2Min(const double &kt2min);
    void SetKT2Max(const double &kt2max);

    void Add(LEV_Base *const lev);

    bool SelectSplitting(const double &rn);

  };// end of class LEV_Group

}// end of namespace EXTRAXS

#endif
