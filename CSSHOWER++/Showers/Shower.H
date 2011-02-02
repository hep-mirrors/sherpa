#ifndef CSSHOWER_Showers_Shower_H
#define CSSHOWER_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Kinematics_Base.H"

namespace ATOOLS { class Cluster_Leg; }

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
    PDF::ISR_Handler * p_isr;
    Parton           * m_last[4];
    int                m_kscheme, m_noem;

    ATOOLS::Cluster_Leg *p_old[2];

    std::map<std::string,double,ATOOLS::String_Sort> m_efac;

    double m_weight;

    void   ResetScales(Parton *const split);
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
    int    RemnantTest(Parton *const p);
    int    UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,int mode);
    int MakeKinematics(Parton *const split,const ATOOLS::Flavour &fla,
		       const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		       const int mode);
 public:
    Shower(PDF::ISR_Handler * isr,const int qed,
	   ATOOLS::Data_Reader *const dataread);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    int    ReconstructDaughters(Singlet *const split,const int mode,
				Parton *const pi=NULL,Parton *const pj=NULL);
    void SetMS(ATOOLS::Mass_Selector *const ms);
    double EFac(const std::string &sfk) const;
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

    inline double Weight() const { return m_weight; }

    inline int KinScheme() const { return m_kscheme; }

  };
}

#endif
