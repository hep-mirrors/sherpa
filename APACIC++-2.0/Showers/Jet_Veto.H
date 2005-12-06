#ifndef Jet_Veto_H
#define Jet_Veto_H

#include "Jet_Finder.H"
#include "Cluster_Algorithm.H"
#include "Timelike_Kinematics.H"
#include "Tree.H"

namespace APACIC {

  class PT_Measure {
  private:
    
    ATOOLS::Jet_Finder *p_jf;

  public:

    double operator()(const ATOOLS::Vec4D &p1);
    double operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2); 

    inline void SetJetFinder(ATOOLS::Jet_Finder *const jf) { p_jf=jf; }

  };// end of class PT_Measure

  class E_Scheme {
  public:
    
    ATOOLS::Vec4D operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2);

  };// end of class E_Scheme

  struct jv {

    enum mode {
      none    = 0,
      final   = 1,
      initial = 2,
      global  = 4,
      mlm     = 8
    };

  };// end of struct jv

  inline jv::mode operator|(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a|(int)b); }
  inline jv::mode operator&(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a&(int)b); }

  class Jet_Veto {
  public:

    typedef ATOOLS::Cluster_Algorithm
    <ATOOLS::Vec4D,PT_Measure,E_Scheme> Cluster_Type;

    typedef std::vector<double> Double_Vector;

  private:

    ATOOLS::Jet_Finder  *p_jf;
    Cluster_Type        *p_cluster;
    Timelike_Kinematics *p_kin;

    Tree **p_istrees, *p_fstree;
    Knot  *p_cur;

    Double_Vector m_rates;

    jv::mode m_mode;

    size_t m_maxjets, m_cmode, m_jmode, m_ljmode;

    double m_ycut;

    int CollectISMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);
    int CollectFSMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);

  public:

    // constructor
    Jet_Veto(ATOOLS::Jet_Finder *const jf,Timelike_Kinematics *const kin);

    //member functions 
    int TestKinematics(const int mode=0,Knot *const mo=NULL);

    int TestFSKinematics(Knot *const mo);
    int TestISKinematics(Knot *const mo);

    // inline functions
    inline void SetYCut(const double &ycut) { m_ycut=ycut; }

    inline void SetMode(const jv::mode &mode)     { m_mode=mode;       }
    inline void SetMaxJets(const size_t &maxjets) { m_maxjets=maxjets; }

    inline void SetJetVeto(const size_t &mode)     { m_jmode=mode;  }
    inline void SetLoseJetVeto(const size_t &mode) { m_ljmode=mode; }

    inline void SetISTrees(Tree **const trees) { p_istrees=trees; }
    inline void SetFSTree(Tree *const tree)    { p_fstree=tree;   }

    inline void SetJetPT2(const double &pt2) { p_jf->SetShowerPt2(pt2); }

    inline jv::mode Mode() const { return m_mode;   }

    inline size_t JetVeto() const     { return m_jmode;  }
    inline size_t LoseJetVeto() const { return m_ljmode; }

    inline ATOOLS::Jet_Finder *const JetFinder() const { return p_jf; }

    inline const Double_Vector &JetRates() const { return m_rates; }

  };// end of class Jet_Veto

}//end of namespace APACIC

#endif
