#ifndef CSSHOWER_Showers_Shower_H
#define CSSHOWER_Showers_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Kinematics_Base.H"

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
    Parton           * m_last[3];

    int    SetXBj(Parton *const p) const;
    PLiter SelectSplitting(double &,const double &);
    bool   EvolveSinglet(Singlet *,const size_t &,size_t &);
    bool   PerformSplitting();
    bool   TrialEmission(double &,const double &,Parton *);
    int    RemnantTest(Parton *const p);
    void   UpdateDaughters(Parton *const split);
    bool   UpdateDaughters(Parton *const split,Parton *const newpB,
			   Parton *const newpC,const bool del=true);
 public:
    Shower(PDF::ISR_Handler * isr,const int qed,
	   ATOOLS::Data_Reader *const dataread);
    ~Shower();
    bool   EvolveShower(Singlet *,const size_t &,size_t &);
    bool   ReconstructDaughters(Singlet *const split,
				const bool one=false);
    void SetMS(ATOOLS::Mass_Selector *const ms);
    void SetRBMax(const double &rbmax);
    inline const Sudakov *GetSudakov() const { return &m_sudakov; }
    inline Parton *const *GetLast() const { return m_last; }

    inline const Kinematics_FF *KinFF() const { return &m_kinFF; }
    inline const Kinematics_FI *KinFI() const { return &m_kinFI; }
    inline const Kinematics_IF *KinIF() const { return &m_kinIF; }
    inline const Kinematics_II *KinII() const { return &m_kinII; }
  };
}

#endif
