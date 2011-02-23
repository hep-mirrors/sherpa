#ifndef PHASIC_Scales_Scale_Setter_Base_H
#define PHASIC_Scales_Scale_Setter_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Vector.H"

namespace MODEL { class Model_Base; }

namespace PHASIC {

  class Process_Base;

  struct stp {

    enum id {
      fac     =  0,
      ren     =  1,
      size    =  2
    };

  };// end of struct stp

  std::ostream &operator<<(std::ostream &ostr,const stp::id &scl);

  class Integrable_Base;
  class Jet_Finder;

  struct Scale_Setter_Arguments {
    Process_Base *p_proc;
    MODEL::Model_Base *p_model;
    MODEL::Coupling_Map *p_cpls;
    std::string m_scale, m_coupling;
    Scale_Setter_Arguments(MODEL::Model_Base *const model,
			   const std::string &scale="",
			   const std::string &coupling=""):
      p_proc(NULL), p_model(model), p_cpls(NULL),
      m_scale(scale), m_coupling(coupling) {}
  };// end of class Scale_Setter_Arguments

  class Scale_Setter_Base {
  public:

    typedef ATOOLS::Getter_Function
    <Scale_Setter_Base,Scale_Setter_Arguments> Scale_Getter_Function;

  protected:

    Process_Base *p_proc, *p_caller;

    MODEL::Model_Base   *p_model;
    MODEL::Coupling_Map *p_cpls;

    std::vector<double> m_scale, m_escale;

    std::string m_coupling;

    bool m_scale2;

    ATOOLS::ClusterAmplitude_Vector m_ampls;

    mutable Jet_Finder *p_jf;

  public:

    Scale_Setter_Base(const Scale_Setter_Arguments &args,
		      const bool scale2=false);

    virtual ~Scale_Setter_Base();

    void SetCouplings();

    virtual double CalculateScale(const std::vector<ATOOLS::Vec4D> &p,
				  const int mode) = 0;

    virtual ATOOLS::Vec4D Momentum(const size_t &i) const;

    virtual double HT() const;

    inline void SetFixedScale(const std::vector<double> &s) { m_escale=s; }

    inline void SetCaller(Process_Base *const proc) { p_caller=proc; }

    inline double Scale(const stp::id type) { return m_scale[type]; }

    inline const std::vector<double> &Scales() const { return m_scale; }

    inline Process_Base *Process() const { return p_proc; }

    inline bool Scale2() { return m_scale2; }

    inline ATOOLS::ClusterAmplitude_Vector &Amplitudes() { return m_ampls; }

    inline MODEL::Coupling_Map *Couplings() const { return p_cpls; }

    static void ShowSyntax(const size_t i);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

#endif
