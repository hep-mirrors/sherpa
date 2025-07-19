#ifndef MALARIC__Shower__Kernel_H
#define MALARIC__Shower__Kernel_H

#include "MALARIC/Shower/Lorentz.H"
#include "MALARIC/Shower/Gauge.H"
#include "MALARIC/Tools/Weight.H"

#include <vector>

namespace MALARIC {

  class Shower;

  class Kernel {
  private:

    Shower *p_ps;

    Lorentz *p_lf;
    Gauge   *p_gf;

    double m_ef;
    int    m_type, m_mode, m_on, m_id;

  public:

    Kernel(Shower *const ps,Kernel_Key key);

    ~Kernel();

    double Value(const Splitting &s) const;

    MC_Weight GetWeight(const Splitting &s,const double &o,
			const MC_Weight *w=NULL) const;

    bool GeneratePoint(Splitting &s) const;

    double Integral(Splitting &s) const;

    int Construct(Splitting &s,const int mode) const;

    std::string Class() const;

    inline bool Allowed(const Splitting &s)
    { return p_lf->Allowed(s)&&p_gf->Allowed(s); }

    inline Lorentz *LF() const { return p_lf; } 
    inline Gauge   *GF() const { return p_gf; } 

    inline Shower *PS() const { return p_ps; }

    inline int Type() const { return m_type; }
    inline int Mode() const { return m_mode; }

    inline int On() const { return m_on; }
    inline int Id() const { return m_id; }

    inline void SetOn(const int &on) { m_on=on; }
    inline void SetId(const int &id) { m_id=id; }

    inline void SetEF(const double &ef) { m_ef=ef; }

  };// end of class Kernel

  typedef std::vector<Kernel*> Kernel_Vector;

}// end of namepsace MALARIC

#endif
