#ifndef METOOLS__Explicit__Dipole_Terms_H
#define METOOLS__Explicit__Dipole_Terms_H

#include <iostream>

namespace METOOLS {

  struct NLO_Value {
    double m_e2, m_e1, m_f;
    inline NLO_Value(const double &e2=0.0,
		     const double &e1=0.0,
		     const double &f=0.0):
      m_e2(e2), m_e1(e1), m_f(f) {}
    inline NLO_Value &operator+=(const NLO_Value &v)
    { m_e2+=v.m_e2; m_e1+=v.m_e1; m_f+=v.m_f; return *this; }
    inline NLO_Value &operator*=(const double &d)
    { m_e2*=d; m_e1*=d; m_f*=d; return *this; }
  };// end of struct NLO_Value

  inline NLO_Value operator+(const NLO_Value &v1,const NLO_Value &v2)
  { return NLO_Value(v1.m_e2+v2.m_e2,v1.m_e1+v2.m_e1,v1.m_f+v2.m_f); }
  inline NLO_Value operator*(const double &d,const NLO_Value &v)
  { return NLO_Value(v.m_e2*d,v.m_e1*d,v.m_f*d); }

  std::ostream &operator<<(std::ostream &str,const NLO_Value &v);

  struct I_Args {
    double Q, Q2, s, mij, mij2, mk, mk2;
    double v, r, r2, rj2, rk2, Qa2;
    I_Args(const double &_s,const double &_mij, const double &_mk);
    void Swap();
  };// end of struct I_Args

  class Dipole_Info;

  NLO_Value FFVS(const I_Args &a,const Dipole_Info *info);

  double FFAE(const I_Args &a,const Dipole_Info *info);

  NLO_Value FFQQ(const I_Args &a,const Dipole_Info *info);
  NLO_Value FFGQQ(const I_Args &a,const Dipole_Info *info);
  double    FFVNSQQ(const I_Args &a,const Dipole_Info *info);
  double    FFACQQ(const I_Args &a,const Dipole_Info *info);

  NLO_Value FFGQ(const I_Args &a,const Dipole_Info *info,const double &m);
  NLO_Value FFGGQ(const I_Args &a,const Dipole_Info *info,const double &m);
  double    FFVNSGQ(const I_Args &a,const Dipole_Info *info,const double &m);
  double    FFACGQ(const I_Args &a,const Dipole_Info *info,const double &m);

  NLO_Value FFGG(const I_Args &a,const Dipole_Info *info);
  NLO_Value FFGGG(const I_Args &a,const Dipole_Info *info);
  double    FFVNSGG(const I_Args &a,const Dipole_Info *info);
  double    FFACGG(const I_Args &a,const Dipole_Info *info);

}// end of namespace METOOLS

#endif
