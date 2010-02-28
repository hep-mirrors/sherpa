#ifndef PHASIC_Scales_Scale_Setter_Base_H
#define PHASIC_Scales_Scale_Setter_Base_H

#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/Vector.H"

namespace PHASIC {

  class Process_Base;

  struct stp {

    enum id {
      ren     =  0,
      fac     =  1,
      size    =  2
    };

  };// end of struct stp

  std::ostream &operator<<(std::ostream &ostr,const stp::id &scl);

  class Integrable_Base;
  class Jet_Finder;

  struct Scale_Setter_Arguments {
    Process_Base *p_proc;
    std::string m_scale;
    Scale_Setter_Arguments(Process_Base *const proc,
			   const std::string &scale=""):
      p_proc(proc), m_scale(scale) {}
  };// end of class Scale_Setter_Arguments

  class Scale_Setter_Base {
  public:

    typedef ATOOLS::Getter_Function
    <Scale_Setter_Base,Scale_Setter_Arguments> Scale_Getter_Function;

  protected:

    Process_Base *p_proc;

    std::vector<double> m_scale;

    ATOOLS::Info_Key m_kfkey;

    bool m_scale2;

    mutable Jet_Finder *p_jf;

  public:

    inline Scale_Setter_Base(Process_Base *const proc,
			     const bool scale2=false):
      p_proc(proc), m_scale(stp::size), m_scale2(scale2), p_jf(NULL) {}

    virtual ~Scale_Setter_Base();

    virtual double CalculateScale(const std::vector<ATOOLS::Vec4D> &p) = 0;
    virtual double CalculateScale2(const std::vector<ATOOLS::Vec4D> &p);

    virtual ATOOLS::Vec4D Momentum(const size_t &i) const;

    virtual double HT() const;
    virtual double YCut() const;

    inline const double Scale(const stp::id type) { return m_scale[type]; }

    inline Process_Base *Process() const { return p_proc; }

    inline bool Scale2() { return m_scale2; }

    static void ShowSyntax(const size_t i);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

#endif
