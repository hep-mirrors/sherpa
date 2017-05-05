#ifndef DIRE__Shower__Kernel_H
#define DIRE__Shower__Kernel_H

#include "DIRE/Shower/Lorentz.H"
#include "DIRE/Shower/Gauge.H"
#include "DIRE/Tools/Weight.H"

#include <vector>

namespace DIRE {

  class Shower;

  class Kernel {
  private:

    Shower *p_ps;

    Lorentz *p_lf;
    Gauge   *p_gf;

    double m_ef;
    int    m_type, m_mode, m_on;

  public:

    Kernel(Shower *const ps,Kernel_Key key);

    ~Kernel();

    double Value(const Splitting &s) const;

    Weight GetWeight(const Splitting &s,const double &o) const;

    bool GeneratePoint(Splitting &s) const;

    double Integral(Splitting &s) const;

    int Construct(Splitting &s,const int mode) const;

    inline bool Allowed(const Splitting &s)
    { return p_lf->Allowed(s)&&p_gf->Allowed(s); }

    inline Lorentz *LF() const { return p_lf; } 
    inline Gauge   *GF() const { return p_gf; } 

    inline Shower *PS() const { return p_ps; }

    inline int Type() const { return m_type; }
    inline int Mode() const { return m_mode; }

    inline int On() const { return m_on; }

    inline void SetEF(const double &ef) { m_ef=ef; }

  };// end of class Kernel

  typedef std::vector<Kernel*> Kernel_Vector;

}// end of namepsace DIRE

#endif
