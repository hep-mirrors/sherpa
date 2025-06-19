#ifndef AMEGIC_DipoleSubtraction_DipoleSplitting_Base_H
#define AMEGIC_DipoleSubtraction_DipoleSplitting_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Recoil_Definition.H"
#include "PDF/Main/NLOMC_Base.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "AMEGIC++/Main/Pol_Info.H"
#include "MODEL/Main/Coupling_Data.H"
#include <vector>
#include <limits>

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

#define nan numeric_limits<double>::quiet_NaN()

  class DipoleSplitting_Base {
  protected:
    std::string m_name;
    double m_alpha,m_kt2max,m_amin,m_kappa;
    double m_Q2,m_kt2,m_a;
    double m_sff,m_av,m_fac;
    int m_mcmode, m_mcsign;
    ATOOLS::Vec4D m_pt1,m_pt2;

    ATOOLS::sbt::subtype       m_stype;
    ATOOLS::dpt::dipoletype    m_dtype;
    ATOOLS::spt::splittingtype m_ftype;
    ATOOLS::Vec4D m_pi,m_pj,m_pk,m_ptij,m_ptk;
    int m_i,m_j,m_k,m_tij,m_tk,m_m,m_es;
    double m_k0sqf,m_k0sqi;
    std::vector<ATOOLS::Vec4D> m_dpollist;
    std::vector<ATOOLS::Vec4D> m_mom;
    std::vector<double> m_pfactors;
    ATOOLS::subscheme::code m_subtype;
    double m_spfdef;
    bool m_collVFF;
    int m_Vsubmode;
    MODEL::Coupling_Data *p_cpl;
    ATOOLS::NLO_subevt *p_subevt;
    PDF::NLOMC_Base *p_nlomc;
    ATOOLS::Recoil_Definition *p_recoil;
    void CalcVectors(ATOOLS::Vec4D&,ATOOLS::Vec4D&);
    void CalcVectors(ATOOLS::Vec4D&,ATOOLS::Vec4D&,double);
  public:
    DipoleSplitting_Base(ATOOLS::sbt::subtype st, ATOOLS::spt::splittingtype ft,
                         ATOOLS::dpt::dipoletype dt,int m,int i,int j,int k);
    virtual ~DipoleSplitting_Base() {}

    inline std::string Name() const                 { return m_name; }

    void SetCoupling(const MODEL::Coupling_Map *cpls);
    inline void   SetChargeFactor(double f) { m_spfdef*=f; }
    inline double SPFac()                   { return m_spfdef*p_cpl->Factor(); }

    virtual void SetMomenta(const ATOOLS::Vec4D*) = 0;
    inline std::vector<ATOOLS::Vec4D>* GetMomenta() { return &m_mom; }
    inline ATOOLS::Vec4D Getptij()                  { return m_ptij; }
    inline ATOOLS::Vec4D Getptk()                   { return m_ptk; }

    inline void SetMCMode(const int mcmode)         { m_mcmode=mcmode; }
    inline void SetAlpha(double alp)                { m_alpha=alp; }
    inline void SetKt2Max(double kt2m)              { m_kt2max=kt2m; }
    inline double Alpha() const                     { return m_alpha; }
    inline double KT2() const                       { return m_kt2; }

    bool Reject(const double &alpha);
    bool KinCheck() const;
    double GetF();
    virtual double GetDPSF() = 0;
    virtual double GetValue() = 0;

    virtual void CalcDiPolarizations() = 0;
    double GetR(const ATOOLS::Vec4D*,const ATOOLS::Vec4D*);
    inline void SetSubevt(ATOOLS::NLO_subevt *const sub) { p_subevt=sub; }
    inline std::vector<ATOOLS::Vec4D>* GetDiPolarizations()
                                                    { return &m_dpollist; }
    inline std::vector<double>* GetFactors()        { return &m_pfactors; }

    inline int MCSign() const                       { return m_mcsign; }
    inline const double& LastAlpha() const          { return m_a; }

    inline const ATOOLS::sbt::subtype       GetSubtractionType() const { return m_stype; }
    inline const ATOOLS::dpt::dipoletype    GetDipoleType() const      { return m_dtype; }
    inline const ATOOLS::spt::splittingtype GetSplittingType() const   { return m_ftype; }
    inline void SetNLOMC(PDF::NLOMC_Base *const mc)
    { p_nlomc=mc; m_subtype=p_nlomc->SubtractionType();
      if (m_subtype==ATOOLS::subscheme::Dire) m_kappa=1.0; }
    inline void SetRecoil(ATOOLS::Recoil_Definition* recoil)
    { p_recoil=recoil; }
  };

  inline std::ostream & operator<<(std::ostream & s,
                                   const DipoleSplitting_Base & dsb)
  {
    return s<<dsb.Name();
  }
}
#endif
