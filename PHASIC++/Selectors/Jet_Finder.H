#ifndef PHASIC_Selectors_Jet_Finder_H
#define PHASIC_Selectors_Jet_Finder_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Phys/Cluster_Leg.H"

namespace PHASIC {
  
  struct Comb_Key: public std::pair<size_t,size_t> {
    size_t partner;
    inline Comb_Key(const size_t &f,const size_t &s,
		    const size_t &p=0):
      std::pair<size_t,size_t>(f,s), partner(p) {}
  };// end of struct Comb_Key

  class Jet_Finder: public Selector_Base {

    double m_ycut, m_gycut, m_cycut, m_gcycut;
    double m_s, m_sprime, m_value;

    std::string m_cuttag, m_gcuttag, m_procname;

    std::map<size_t,ATOOLS::Vec4D>   m_moms;
    std::map<size_t,ATOOLS::Flavour> m_flavs;
    std::map<size_t,ATOOLS::ColorID> m_cols;

    ATOOLS::Flavour_Vector           m_strongflavs;
    std::vector<size_t>              m_stronglocs;

    std::map<size_t,std::map<size_t,double> > m_ycuts, m_gycuts;

    std::vector<std::vector<int> > m_mcomb;

    std::vector<std::vector<Comb_Key> > m_fills;

    bool m_sok, m_on;

    ATOOLS::Flavour GetFlavour(std::string fl);

    size_t FillCombinations(const std::string &name,const std::string &ycut,
			    const std::string &gycut,size_t &cp,const int fl);

    void PrepareMomList(const ATOOLS::Vec4D_Vector &vec);
    bool PrepareColList(const std::vector<int> &ci,
			const std::vector<int> &cj);

    bool ColorConnected(const size_t &i,const size_t &j) const;

  public:

    Jet_Finder(const int nin=0,const int nout=0,ATOOLS::Flavour *const fl=NULL,
	       const std::string &ycut="1",const std::string &gycut="1");

    ~Jet_Finder();

    void   FillCombinations();
    bool   Trigger(const ATOOLS::Vec4D_Vector &p);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,
                      const ATOOLS::Flavour_Vector &, int);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &p);
    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &p);

    double Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
		const ATOOLS::Vec4D &pk,
		const ATOOLS::Flavour &fi=ATOOLS::Flavour(kf_none),
		const ATOOLS::Flavour &fj=ATOOLS::Flavour(kf_none));

    void   BuildCuts(Cut_Data *cuts);
    void   UpdateCuts(double,double,Cut_Data *cuts);

    double ActualValue() const;

    double GetYcut(const size_t &i,const size_t &j) const;
    double GetGlobalYcut(const size_t &i,const size_t &j) const;

    inline double Ycut() const           { return m_ycut;   }
    inline double GlobalYcut() const     { return m_gycut;  }
    inline double CoreYcut() const       { return m_cycut;  }
    inline double GlobalCoreYcut() const { return m_gcycut; }

    inline void SetOn(const bool on) { m_on=on; }

  };

}

#endif








