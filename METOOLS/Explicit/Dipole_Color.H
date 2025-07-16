#ifndef METOOLS_Explicit_Dipole_Color_H
#define METOOLS_Explicit_Dipole_Color_H

#include "METOOLS/Explicit/Color_Calculator.H"

namespace METOOLS {

  class Dipole_Color: public Color_Calculator {
  public:

    struct CInfo {
      int m_cr, m_ca, m_i, m_t;
      double m_s;
      inline CInfo(const int cr,const int ca,
		   const int i,const int t,
		   const double &s=1.0):
	m_cr(cr), m_ca(ca), m_i(i), m_t(t), m_s(s) {}
    };// end of struct CInfo

  protected:

    Dipole_Info *p_dinfo;
    Current     *p_kt;

    Color_Calculator *p_cc;

    std::vector<CInfo> m_ci, m_cjk;

  public:

    // constructor
    Dipole_Color(const Vertex_Key &key);

    // destructor
    ~Dipole_Color();

    // member functions
    void AddJ(CObject *const j);

    void AddJI(CObject *const j,const int t) const;
    void AddJJK(CObject *const j) const;

  };// end of class Dipole_Color

  std::ostream &operator<<(std::ostream &str,
			   const Dipole_Color::CInfo &c);

}// end of namespace METOOLS

#endif
