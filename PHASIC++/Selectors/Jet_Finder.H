#ifndef PHASIC_Selectors_Jet_Finder_H
#define PHASIC_Selectors_Jet_Finder_H

#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Phys/Cluster_Leg.H"

namespace PHASIC {

  class Single_Process;
  
  struct Comb_Key: public std::pair<size_t,size_t> {
    size_t partner;
    int cpl;
    inline Comb_Key(const int &_cpl, const size_t &f,const size_t &s,
		    const size_t &p=0):
      std::pair<size_t,size_t>(f,s), partner(p), cpl(_cpl) {}
  };// end of struct Comb_Key

  class Jet_Finder: public Selector_Base,
		    public ATOOLS::Tag_Replacer {

    double m_ycut, m_s, m_pt2min, m_dparam;

    std::string m_cuttag;

    std::map<size_t,ATOOLS::Vec4D>   m_moms;
    std::map<size_t,ATOOLS::Flavour> m_flavs;
    std::map<size_t,ATOOLS::ColorID> m_cols;

    ATOOLS::Flavour_Vector           m_strongflavs;
    std::vector<size_t>              m_stronglocs;

    std::vector<std::vector<int> > m_mcomb;

    std::vector<std::vector<Comb_Key> > m_fills;

    std::set<size_t> m_pcs;

    bool m_on;

    Single_Process *p_sproc;

    ATOOLS::Algebra_Interpreter *p_yccalc;

    ATOOLS::Flavour GetFlavour(std::string fl);

    size_t FillCombinations(const Subprocess_Info &pinfo,
                            size_t &cp,const int fl);

    void PrepareMomList(const ATOOLS::Vec4D_Vector &vec);
    bool PrepareColList(const std::vector<int> &ci,
			const std::vector<int> &cj);

    bool ColorConnected(const size_t &i,const size_t &j,const size_t &k) const;
    bool ColorConnected(const ATOOLS::ColorID &ci,const ATOOLS::ColorID &cj,
			const int si,const int sj) const;

  public:

    Jet_Finder(const int nin=0,const int nout=0,ATOOLS::Flavour *const fl=NULL,
	       const std::string &ycut="1");

    ~Jet_Finder();

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
    void          AssignId(ATOOLS::Term *term);

    void   FillCombinations();
    bool   Trigger(const ATOOLS::Vec4D_Vector &p);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,
                      const ATOOLS::Flavour_Vector &, int);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &p);
    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &p);

    double Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
		const ATOOLS::Vec4D &pk,
		const ATOOLS::Flavour &fi=ATOOLS::Flavour(kf_none),
		const ATOOLS::Flavour &fj=ATOOLS::Flavour(kf_none),
		const int mode=0);

    void   BuildCuts(Cut_Data *cuts);
    void   UpdateCuts(double,double,Cut_Data *cuts);

    inline double Ycut() const           { return m_ycut;   }

    inline void SetOn(const bool on) { m_on=on; }

  };

}

#endif








