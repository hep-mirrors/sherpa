#ifndef PHASIC_Process_Virtual_ME2_Base_H
#define PHASIC_Process_Virtual_ME2_Base_H

#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "MODEL/Main/Model_Base.H"

namespace PHASIC {

  class Process_Info;

  class Virtual_ME2_Base {
  protected:
    std::string m_name;
    const Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    ATOOLS::sbt::subtype m_stype;
    METOOLS::DivArrD m_res;
    double m_born, m_norm;
    double m_mur2, m_accu;
    int m_mode, m_drmode, m_colmode;
    bool m_providespoles, m_fixedIRscale;
    double m_IRscale, m_UVscale;
    MODEL::Coupling_Data* p_aqcd, * p_aqed;
    std::vector<double> m_asscontribs;
  public:
    Virtual_ME2_Base(const PHASIC::Process_Info& pi,
                  const ATOOLS::Flavour_Vector& flavs);
    virtual ~Virtual_ME2_Base();

    virtual void SwitchMode(const int mode=0);

    // Allow to pass a born value for virtuals that don't calculate it
    // themselves but need it for the renormalization scale dependence
    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta,
		      const double& born2) { return Calc(momenta); }

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta)=0;
    
    inline const METOOLS::DivArrD& Result() { return m_res; }

    virtual double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    virtual double ScaleDependenceCoefficient(const int i);
    virtual bool   IsMappableTo(const Process_Info& pi);

    inline double ME_Finite() { return m_res.GetFinite(); }
    inline double ME_E1() { return m_res.GetIR(); }
    inline double ME_E2() { return m_res.GetIR2(); }
    inline double ME_Born() { return m_born; }
    
    inline bool ProvidesPoles() const { return m_providespoles; }
    inline bool fixedIRscale() const  { return m_fixedIRscale; }
    inline double IRscale() const { return m_IRscale; }
    inline double UVscale() const { return m_UVscale; }

    inline size_t ME_AssContribs_Size() { return m_asscontribs.size(); }
    inline double ME_AssContribs(size_t i) { return m_asscontribs[i]; }

    inline void SetSubType(const ATOOLS::sbt::subtype& st) { m_stype=st; }
    inline void SetRenScale(const double& mur2) { m_mur2=mur2; }
    inline void SetNorm(const double& norm) { m_norm=norm; }
    void SetCouplings(const MODEL::Coupling_Map& cpls);

    double AlphaQCD() const;
    double AlphaQED() const;

    inline int Mode() const    { return m_mode;    }
    inline int DRMode() const  { return m_drmode;  }
    inline int ColMode() const { return m_colmode; }

    inline std::string Name() { return m_name; }

    static Virtual_ME2_Base* GetME2(const PHASIC::Process_Info& pi);
    static Virtual_ME2_Base* GetME2(const std::string& tag,
                                    const PHASIC::Process_Info& pi);
  };

  class Trivial_Virtual : public Virtual_ME2_Base {
  public:
    Trivial_Virtual(const PHASIC::Process_Info& pi,
                 const ATOOLS::Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs) {}

    inline void Calc(const ATOOLS::Vec4D_Vector& momenta) {}
  };
}

#define DECLARE_VIRTUALME2_GETTER(NAME,TAG)			            \
  DECLARE_GETTER(NAME,TAG,PHASIC::Virtual_ME2_Base,PHASIC::Process_Info);   \
  void ATOOLS::Getter<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,NAME>:: \
  PrintInfo(std::ostream &str,const size_t width) const		            \
  {                                                                         \
    str<<#TAG;                                                              \
  }

#endif
