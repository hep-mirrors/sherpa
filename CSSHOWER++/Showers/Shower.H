#ifndef CSSHOWER_Showers_Shower_H
#define CSSHOWER_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Phys/Weights.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Kinematics_Base.H"

namespace ATOOLS {
  class Cluster_Leg;
}

namespace CSSHOWER {

  class Shower {
  private:
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    Kinematics_II      m_kinII;
    ATOOLS::Flavour    m_flavA, m_flavB, m_flavC;
    Singlet      *     p_actual;
    Sudakov            m_sudakov;
    PDF::ISR_Handler          * p_isr;
    REMNANTS::Remnant_Handler * p_remnants;
    Parton           * m_last[4];
    int                m_kscheme, m_maxpart, m_use_bbw, m_recdec;
    bool               m_mecorrections;

    ATOOLS::Cluster_Leg *p_old[2];

    std::map<std::string,double,ATOOLS::String_Sort> m_efac;

    ATOOLS::Weights_Map m_weightsmap;

    bool m_reweight;
    double m_maxreweightfactor;

    double m_lastcpl;

    void   ResetScales(const double &kt2);
    void   SetSplitInfo(const ATOOLS::Vec4D &psplit,
			const ATOOLS::Vec4D &pspect,
			Parton *const split,Parton *const newb,
			Parton *const newc,const int mode);
    double GetXBj(Parton *const p) const;
    int    SetXBj(Parton *const p) const;
    Parton *SelectSplitting(double &);
    void   EnableMECorrection(Singlet *actual);
    bool   EvolveSinglet(Singlet *,const size_t &,size_t &);
    bool   PerformSplitting();
    bool   TrialEmission(double &,Parton *);
    int    RemnantTest(Parton *const p,const ATOOLS::Poincare_Sequence *lt);
    int    UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,double &jcv,int mode);
    int MakeKinematics(Parton *const split,const ATOOLS::Flavour &fla,
		       const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		       double &jcv,int mode);

  public:
    Shower(PDF::ISR_Handler *,const int qcd,const int qed, int type);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    int    ReconstructDaughters(Singlet *const split,double &jcv,
				Parton *const pi,Parton *const pj);
    double Reweight(ATOOLS::QCD_Variation_Params*, Parton& splitter);
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

    inline const ATOOLS::Weights_Map& WeightsMap() const { return m_weightsmap; }

    inline int KinScheme() const { return m_kscheme; }

    inline PDF::ISR_Handler *ISR() { return p_isr; }

    inline void SetRemnants(REMNANTS::Remnant_Handler * remnants) {
      p_remnants = remnants;
    }
  };
}

#endif
