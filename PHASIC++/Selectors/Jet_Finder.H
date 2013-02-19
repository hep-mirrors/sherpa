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

    double m_ycut, m_s, m_dparam;
    std::string m_cuttag;
    bool m_on;

    ATOOLS::Algebra_Interpreter *p_yccalc;

    ATOOLS::Cluster_Amplitude *p_ampl;

    PDF::Jet_Criterion *p_jc;

    ATOOLS::Flavour GetFlavour(std::string fl);

  public:

    Jet_Finder(const int nin=0,const int nout=0,
	       ATOOLS::Flavour *const fl=NULL,
	       const std::string &ycut="1");

    ~Jet_Finder();

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
    void          AssignId(ATOOLS::Term *term);

    bool   Trigger(const ATOOLS::Vec4D_Vector &p);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &p,
		      ATOOLS::NLO_subevtlist *const subs);
    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &p);

    void   BuildCuts(Cut_Data *cuts);

    inline double Ycut() const           { return m_ycut;   }
    inline double DR() const             { return m_dparam; }

    inline void SetOn(const bool on) { m_on=on; }

    inline PDF::Jet_Criterion *JC() const { return p_jc; }

  };

}

#endif








