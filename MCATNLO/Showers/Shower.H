#ifndef MCATNLO_Showers_Shower_H
#define MCATNLO_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Phys/Event_Weights.H"
#include "MCATNLO/Showers/Sudakov.H"
#include "MCATNLO/Tools/Singlet.H"
#include "MCATNLO/Showers/Kinematics_Base.H"

namespace ATOOLS {

  class Cluster_Leg;

}

namespace MCATNLO {

  class CS_Gamma;

  class Shower {
  private:
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    Kinematics_II      m_kinII;
    ATOOLS::Flavour    m_flavA, m_flavB, m_flavC, m_flav;
    Singlet      *     p_actual;
    Sudakov            m_sudakov;
    PDF::ISR_Handler * p_isr;
    Parton           * m_last[4];
    int                m_kscheme;

    ATOOLS::Cluster_Leg *p_old[2];

    CS_Gamma *p_gamma;
    ATOOLS::Event_Weights m_weights;

    /// whether to disable the reweighting of the emission
    bool m_norewem;

    void   ResetScales(Parton *const split);
    void   SetSplitInfo(const ATOOLS::Vec4D &psplit,
			const ATOOLS::Vec4D &pspect,
			Parton *const split,Parton *const newb,
			Parton *const newc,const int mode);
    int    SetXBj(Parton *const p) const;
    Parton *SelectSplitting(double &);
    bool   PerformSplitting();
    bool   TrialEmission(double &,Parton *);
    int    RemnantTest(Parton *const p);
    int    UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,double &jcv);
    int MakeKinematics(Parton *const split,const ATOOLS::Flavour &fla,
		       const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		       double &jcv);

  public:
    Shower(PDF::ISR_Handler *,const int qed);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    void SetMS(const ATOOLS::Mass_Selector *const ms);
    inline Sudakov *GetSudakov() { return &m_sudakov; }
    inline const Sudakov *GetSudakov() const { return &m_sudakov; }
    inline Parton *const *GetLast() const { return m_last; }
    inline ATOOLS::Cluster_Leg *GetOld(const int i) const { return p_old[i]; }

    inline const Kinematics_FF *KinFF() const { return &m_kinFF; }
    inline const Kinematics_FI *KinFI() const { return &m_kinFI; }
    inline const Kinematics_IF *KinIF() const { return &m_kinIF; }
    inline const Kinematics_II *KinII() const { return &m_kinII; }

    inline Kinematics_FF *KinFF() { return &m_kinFF; }
    inline Kinematics_FI *KinFI() { return &m_kinFI; }
    inline Kinematics_IF *KinIF() { return &m_kinIF; }
    inline Kinematics_II *KinII() { return &m_kinII; }

    inline int KinScheme() const { return m_kscheme; }

    inline void SetGamma(CS_Gamma *const gamma) { p_gamma=gamma; }

    inline CS_Gamma *Gamma() const { return p_gamma; }

    inline ATOOLS::Event_Weights Weights() const { return m_weights; }

    inline const ATOOLS::Flavour &ActiveFlav() const { return m_flav; }

    inline PDF::ISR_Handler *ISR() { return p_isr; }
  };
}

#endif
