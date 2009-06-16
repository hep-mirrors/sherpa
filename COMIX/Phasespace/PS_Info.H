#ifndef COMIX_Phasespace_PS_Info_H
#define COMIX_Phasespace_PS_Info_H

#include <iostream>

namespace COMIX {

  class PS_Info {
  private:

    double m_w;

    int m_c[2];

    long unsigned int m_h;

    friend std::ostream &operator<<(std::ostream &str,const PS_Info &s);

  public:

    inline PS_Info(const PS_Info &i): m_w(i.m_w), m_h(i.m_h)
    { 
      m_c[0]=i.m_c[0]; m_c[1]=i.m_c[1];
    }
    inline PS_Info(const int &c1=0,const int &c2=0,
		   const long unsigned int &h=0,
		   const double &w=1.0): m_w(w), m_h(h)
    { 
      m_c[0]=c1; m_c[1]=c2;
    }

    inline int    &operator()(const int i) { return m_c[i]; }
    inline double &operator[](const int i) { return m_w;    }

    inline const int    &operator()(const int i) const { return m_c[i]; }
    inline const double &operator[](const int i) const { return m_w;    }

    inline void SetH(const long unsigned int &h) { m_h=h; }

    inline long unsigned int H() const { return m_h; }

    inline PS_Info operator-() const { return PS_Info(m_c[0],m_c[1],m_h); }

    inline bool operator==(const PS_Info &i) const
    {
      return m_h==i.m_h && m_c[0]==i.m_c[0] && m_c[1]==i.m_c[1];
    }

  };// end of class PS_Info

  std::ostream &operator<<(std::ostream &str,const PS_Info &s);

}// end of namespace ATOOLS

#endif
