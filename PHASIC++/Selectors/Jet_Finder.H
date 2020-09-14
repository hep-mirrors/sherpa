#ifndef PHASIC_Selectors_Jet_Finder_H
#define PHASIC_Selectors_Jet_Finder_H

#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PDF { class Jet_Criterion; }

namespace PHASIC {

  class Jet_Finder: public Selector_Base,
		    public ATOOLS::Tag_Replacer {

    struct ReweightSubevt_Args {
      std::vector<double> m_jcv;
      std::vector<int> m_acc;
      ReweightSubevt_Args(const size_t &n):
	m_jcv(n,0.0), m_acc(n,0) {}
    };// end of struct ReweightSubevt_Args

    double m_qcut;
    std::string m_cuttag;
    bool m_on;

    ATOOLS::Algebra_Interpreter *p_yccalc;
    ATOOLS::Cluster_Amplitude *p_ampl;
    ATOOLS::Flavour GetFlavour(std::string fl);

  protected:
    PDF::Jet_Criterion *p_jc;
  public:
    Jet_Finder() : Selector_Base("empty_JF"), m_qcut(1.),
		   p_yccalc(NULL), p_ampl(NULL), p_jc(NULL) {};
    Jet_Finder(Process_Base *const proc,const std::string &ycut="1");

    ~Jet_Finder();

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
    void          AssignId(ATOOLS::Term *term);

    bool   Trigger(ATOOLS::Selector_List &p);
    bool   RSTrigger(ATOOLS::NLO_subevtlist *const subs);

    void   BuildCuts(Cut_Data *cuts);

    inline double Qcut() const { return m_qcut; }

    inline void SetOn(const bool on) { m_on=on; }

    inline PDF::Jet_Criterion *JC() const { return p_jc; }

  };

}

#endif








