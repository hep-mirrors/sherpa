#ifndef Jet_Veto_H
#define Jet_Veto_H

#include "Jet_Finder.H"
#include "Cluster_Algorithm.H"
#include "Timelike_Kinematics.H"
#include "Histogram.H"
#include "Exception.H"
#include "Tree.H"
/*!
  \file Jet_Veto.H
  \brief Declares the class Jet_Veto
*/

#define JET_MEASURE PT_Measure
#define RECOMBINATION_SCHEME P_Scheme

namespace APACIC {

  class PT_Measure {
  private:
    
    ATOOLS::Jet_Finder *p_jf;

  public:

    double operator()(const ATOOLS::Vec4D &p1);
    double operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2); 

    inline void SetJetFinder(ATOOLS::Jet_Finder *const jf) { p_jf=jf; }

  };// end of class PT_Measure
  /*!
    \class PT_Measure
    \brief Distance measure for the Cluster_Algorithm of the Jet_Veto

    Defines the distance between two momenta according to the Jet_Finder.
  */

  class E_Scheme {
  public:
    
    ATOOLS::Vec4D operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2);

    void operator()(const ATOOLS::Vec4D &p1);

  };// end of class E_Scheme
  /*!
    \class E_Scheme
    \brief Recombination scheme for the Cluster_Algorithm of the Jet_Veto

    Combines two momenta \f$p_i\f$ and \f$p_j\f$ by 4-vector addition.
  */

  class P_Scheme {
  public:
    
    ATOOLS::Vec4D operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2);

    void operator()(const ATOOLS::Vec4D &p1);

  };// end of class P_Scheme
  /*!
    \class P_Scheme
    \brief Recombination scheme for the Cluster_Algorithm of the Jet_Veto

    Combines two momenta \f$p_i\f$ and \f$p_j\f$ via
    \f[
    \vec{p}_{ij} = \vec{p}_i+\vec{p}_j\;,
    \f]
    \f[
    E_{ij} = |\vec{p}_i+\vec{p}_j|\;.
    \f]
  */

  struct jv {

    enum mode {
      none    = 0,
      final   = 1,
      initial = 2,
      global  = 4,
      mlm     = 8
    };

  };// end of struct jv
  /*!
    \struct jv
    \brief Declares abbreviations for jet veto schemes.
    
    There are essentially three different methods to do the jet veto.
    The first one is a local jet veto, specified by the final and 
    initial flag, but without the global flag. It means taking into 
    account only the distance between two daughter partons, when 
    testing an emission. The global jet veto, specified by the initial 
    and final and the global flag means redetermining the distances 
    between all partons and rejecting an emission, if the jet number 
    is modified by the current emission. The mlm flag means that the
    jet configuration is tested only once, after the whole shower 
    evolution and the event is rejected as a whole, if the jet number 
    is modified.
  */

  inline jv::mode operator|(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a|(int)b); }
  inline jv::mode operator&(const jv::mode &a,const jv::mode &b)
  { return (jv::mode)((int)a&(int)b); }

  class Jet_Veto: public ATOOLS::Terminator_Object {
  public:

    typedef ATOOLS::Cluster_Algorithm
    <ATOOLS::Vec4D,JET_MEASURE,RECOMBINATION_SCHEME> Cluster_Type;

    typedef std::vector<double> Double_Vector;

    typedef std::vector<ATOOLS::Histogram*> Histogram_Vector;

  private:

    ATOOLS::Jet_Finder  *p_jf;
    Cluster_Type        *p_cluster;
    Timelike_Kinematics *p_kin;

    Tree **p_istrees, *p_fstree;
    Knot  *p_cur;

    Double_Vector m_rates;

    jv::mode m_mode;

    size_t m_maxjets, m_cmode, m_jmode, m_ljmode;

    double m_ycut, m_q2hard;

    Histogram_Vector m_histos;

    int CollectISMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);
    int CollectFSMomenta(Knot *knot,std::vector<ATOOLS::Vec4D> &vecs,
			 size_t &hard);

  public:

    // constructor
    Jet_Veto(ATOOLS::Jet_Finder *const jf,Timelike_Kinematics *const kin);

    // destructor
    ~Jet_Veto();

    //member functions 
    void PrepareTerminate();

    int TestKinematics(const int mode=0,Knot *const mo=NULL);
    /*!
      \fn int TestKinematics(const int mode=0,Knot *const mo=NULL)
      \brief Implements the global jet veto procedure

      This procedure is called, to veto parton configurations according 
      to the number of jets, which they define. First the current final state 
      partons are extracted from p_istrees and p_fstree, then p_cluster is 
      applied to define jets and finally the configuration is rejected 
      according to the settings in m_mode.
    */

    int TestFSKinematics(Knot *const mo);
    /*!
      \fn int TestFSKinematics(Knot *const mo)
      \brief Implements the jet veto for final state splittings.
      
      This method implements the jet veto for final state splittings.
      According to the settings in m_mode, either a local or a global 
      jet veto is applied.
    */
    int TestISKinematics(Knot *const mo);
    /*!
      \fn int TestISKinematics(Knot *const mo)
      \brief Implements the jet veto for initial state splittings.
      
      This method implements the jet veto for initial state splittings.
      According to the settings in m_mode, either a local or a global 
      jet veto is applied.
    */

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
  /*!
    \class Jet_Veto
    \brief Implements the jet veto procedure.

    This class implements the jet veto procedure. It is called by the
    Final_State_Shower, any time a new emission is constructed with the 
    full kinematics and tests the jet measure of this emission w.r.t.
    other partons. The emission is rejected, if a new jet has been created. 
    Likewise, the parton configuration is vetoed, if the number of jets
    drops below the number of initial jets, i.e. this which has been set by 
    the hard matrix element. This is a final check, which is performed 
    only once, after the complete shower evolution. Jet_Veto is then 
    called by Apacic.
  */

}//end of namespace APACIC

#endif
