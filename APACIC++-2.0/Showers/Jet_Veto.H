#ifndef Jet_Veto_H
#define Jet_Veto_H

#include "Jet_Finder.H"
#include "Timelike_Kinematics.H"
#include "Tree.H"

namespace APACIC {

  struct jv {

    enum mode {
      none    = 0,
      final   = 1,
      initial = 2,
      mlm     = 4
    };

  };// end of struct jv

  inline jv::mode operator|(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a|(int)b); }
  inline jv::mode operator&(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a&(int)b); }

  class Jet_Veto {
  public:

    typedef std::vector<double> Double_Vector;

  private:

    ATOOLS::Jet_Finder  *p_jf;
    Timelike_Kinematics *p_kin;

    Tree **p_istrees, *p_fstree;

    Double_Vector m_rates;

    jv::mode m_mode;

    size_t m_maxjets;

    int CollectISMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);
    int CollectFSMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);

  public:

    // constructor
    Jet_Veto(ATOOLS::Jet_Finder *const jf,Timelike_Kinematics *const kin);

    //member functions 
    int  TestKinematics(const int mode=0,Knot *const mo=NULL);

    int TestFSKinematics(Knot *const mo);
    int TestISKinematics(Knot *const mo);

    // inline functions
    inline void SetMode(const jv::mode &mode)     { m_mode=mode;       }
    inline void SetMaxJets(const size_t &maxjets) { m_maxjets=maxjets; }
    
    inline void SetISTrees(Tree **const trees) { p_istrees=trees; }
    inline void SetFSTree(Tree *const tree)    { p_fstree=tree;   }

    inline void SetJetPT2(const double &pt2) { p_jf->SetShowerPt2(pt2); }

    inline jv::mode Mode() const { return m_mode; }

    inline ATOOLS::Jet_Finder *const JetFinder() const { return p_jf; }

    inline const Double_Vector &JetRates() const { return m_rates; }

  };// end of class Jet_Veto

}//end of namespace APACIC

#endif
