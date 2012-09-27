#ifndef PHASIC_Scales_Scale_Setter_Base_H
#define PHASIC_Scales_Scale_Setter_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Vector.H"

namespace MODEL { class Model_Base; }

namespace PHASIC {

  class Process_Base;

  struct Scale_Setter_Arguments {
    Process_Base *p_proc;
    MODEL::Model_Base *p_model;
    MODEL::Coupling_Map *p_cpls;
    std::string m_scale, m_coupling;
    size_t m_nin, m_nout;
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

    ATOOLS::Vec4D_Vector    m_p;
    ATOOLS::NLO_subevtlist *p_subs;

    std::vector<double> m_scale, m_escale, m_fac;

    std::string m_coupling;

    size_t m_nin, m_nout;

    double m_htyfac, m_htyexp;

    ATOOLS::ClusterAmplitude_Vector m_ampls;

    virtual double Calculate(const ATOOLS::Vec4D_Vector &p,
			     const size_t &mode) = 0;

  public:

    Scale_Setter_Base(const Scale_Setter_Arguments &args);

    virtual ~Scale_Setter_Base();

    void SetCouplings();

    double HT() const;

    ATOOLS::Vec4D PSum() const;

    double HTYweighted() const;
    bool   SetHTYweightedParameters(std::string par, std::string val);
    
    double BeamThrust() const;

    double CalculateScale(const ATOOLS::Vec4D_Vector &p,const size_t mode=0);

    inline void SetFixedScale(const std::vector<double> &s) { m_escale=s; }

    inline void SetCaller(Process_Base *const proc) { p_caller=proc; }

    inline double Scale(const ATOOLS::stp::id type,const int mode=0) const
    { return mode==0?m_scale[type]:m_fac[type]*m_scale[type]; }

    inline const std::vector<double> &Scales() const { return m_scale; }

    inline const std::vector<double> &FixedScales() const { return m_escale; }

    inline Process_Base *Process() const { return p_proc; }

    inline ATOOLS::ClusterAmplitude_Vector &Amplitudes() { return m_ampls; }

    inline MODEL::Coupling_Map *Couplings() const { return p_cpls; }

    static void ShowSyntax(const size_t i);

    inline const ATOOLS::Vec4D_Vector &Momenta() const { return m_p; }

    inline size_t NIn() const  { return m_nin;  }
    inline size_t NOut() const { return m_nout; }

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

#endif
