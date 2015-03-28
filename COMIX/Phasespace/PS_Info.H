#ifndef COMIX_Phasespace_PS_Info_H
#define COMIX_Phasespace_PS_Info_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <vector>
#include <iostream>

namespace COMIX {

  class PS_Info: public METOOLS::CObject {
  private:

    double m_w;

    static ATOOLS::AutoDelete_Vector<PS_Info> s_objects;

    friend std::ostream &operator<<(std::ostream &str,const PS_Info &s);

  public:

    static PS_Info *New();
    static PS_Info *New(const PS_Info &s);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    inline PS_Info(const PS_Info &i): m_w(i.m_w)
    { 
      m_c[0]=i.m_c[0]; m_c[1]=i.m_c[1]; m_h=i.m_h; m_s=i.m_s;
    }
    inline PS_Info(const int &c1=0,const int &c2=0,
		   const double &w=1.0): m_w(w)
    { 
      m_c[0]=c1; m_c[1]=c2; m_h=0; m_s=0;
    }

    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Conjugate();
    void Invert();

    inline double &operator[](const int i) { return m_w; }

    inline const double &operator[](const int i) const { return m_w; }

  };// end of class PS_Info

  std::ostream &operator<<(std::ostream &str,const PS_Info &s);

}// end of namespace ATOOLS

#endif
