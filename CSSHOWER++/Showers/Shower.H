#ifndef CSSHOWER_Showers_Shower_H
#define CSSHOWER_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Kinematics_Base.H"

namespace ATOOLS {

  class Cluster_Leg;
  class Variation_Parameters;
  class Variation_Weights; 

}

namespace CSSHOWER {

  struct JetVeto_Args {
    Singlet *p_sing;
    double m_jcv;
    int m_acc, m_nlo;
    std::vector<int> m_skip;
    JetVeto_Args(Singlet *const sing,
		 const double &jcv,const size_t &n):
      p_sing(sing), m_jcv(jcv), m_acc(0), m_skip(n,0) {}
  };// end of struct JetVeto_Args
    
  class Shower {
  private:
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    Kinematics_II      m_kinII;
    ATOOLS::Flavour    m_flavA, m_flavB, m_flavC;
    Singlet      *     p_actual;
    Sudakov            m_sudakov;
    PDF::ISR_Handler * p_isr;
    Parton           * m_last[4];
    int                m_kscheme, m_maxpart, m_use_bbw, m_recdec;

    ATOOLS::Cluster_Leg *p_old[2];

    std::map<std::string,double,ATOOLS::String_Sort> m_efac;

    double m_weight, m_lastcpl;

    ATOOLS::Variation_Weights * p_variationweights;
    //! the maximum number of rew'd emissions
    int m_maxrewem;

    void   ResetScales(const double &kt2);
    void   SetSplitInfo(const ATOOLS::Vec4D &psplit,
			const ATOOLS::Vec4D &pspect,
			Parton *const split,Parton *const newb,
			Parton *const newc,const int mode);
    double GetXBj(Parton *const p) const;
    int    SetXBj(Parton *const p) const;
    Parton *SelectSplitting(double &);
    bool   EvolveSinglet(Singlet *,const size_t &,size_t &);
    bool   PerformSplitting();
    bool   TrialEmission(double &,Parton *);
    int    RemnantTest(Parton *const p,const ATOOLS::Poincare_Sequence *lt);
    int    UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,double &jcv,int mode);
    int MakeKinematics(Parton *const split,const ATOOLS::Flavour &fla,
		       const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		       double &jcv,int mode);
    double GetWeight(ATOOLS::Variation_Parameters *params,
		     ATOOLS::Variation_Weights *weights,
		     std::vector<double> &v);
    double VetoWeight(ATOOLS::Variation_Parameters *params,
		      ATOOLS::Variation_Weights *weights,
		      JetVeto_Args &a);
 public:
    Shower(PDF::ISR_Handler *,const int qcd,const int qed,
	   ATOOLS::Default_Reader *const,int type);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    int    ReconstructDaughters(Singlet *const split,double &jcv,
				Parton *const pi,Parton *const pj);
    void SetMS(const ATOOLS::Mass_Selector *const ms);
    double EFac(const std::string &sfk) const;
    inline Sudakov *GetSudakov() { return &m_sudakov; }
    inline const Sudakov *GetSudakov() const { return &m_sudakov; }

    inline const Kinematics_FF *KinFF() const { return &m_kinFF; }
    inline const Kinematics_FI *KinFI() const { return &m_kinFI; }
    inline const Kinematics_IF *KinIF() const { return &m_kinIF; }
    inline const Kinematics_II *KinII() const { return &m_kinII; }

    inline Kinematics_FF *KinFF() { return &m_kinFF; }
    inline Kinematics_FI *KinFI() { return &m_kinFI; }
    inline Kinematics_IF *KinIF() { return &m_kinIF; }
    inline Kinematics_II *KinII() { return &m_kinII; }

    inline double Weight() const { return m_weight; }

    inline int KinScheme() const { return m_kscheme; }

    inline PDF::ISR_Handler *ISR() { return p_isr; }

    inline void SetVariationWeights(ATOOLS::Variation_Weights * varwgts)
    { p_variationweights = varwgts; }

  };
}

#endif
