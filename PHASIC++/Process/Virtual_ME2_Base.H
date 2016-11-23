#ifndef PHASIC_Process_Virtual_ME2_Base_H
#define PHASIC_Process_Virtual_ME2_Base_H

#include "PHASIC++/Process/Process_Info.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Model_Base.H"

namespace PHASIC {

  class Virtual_ME2_Base {
  protected:
    const Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    METOOLS::DivArrD m_res;
    double m_born, m_norm;
    double m_mur2;
    int m_mode, m_drmode, m_colmode;
    MODEL::Coupling_Data* p_aqcd, * p_aqed;
  public:
    Virtual_ME2_Base(const PHASIC::Process_Info& pi,
                  const ATOOLS::Flavour_Vector& flavs);
    virtual ~Virtual_ME2_Base();

    virtual void SwitchMode(const int mode=0);

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta)=0;
    inline const METOOLS::DivArrD& Result() { return m_res; }

    virtual double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    virtual double ScaleDependenceCoefficient(const int i);

    inline double ME_Finite() { return m_res.GetFinite(); }
    inline double ME_E1() { return m_res.GetIR(); }
    inline double ME_E2() { return m_res.GetIR2(); }
    inline double ME_Born() { return m_born; }

    inline void SetRenScale(const double& mur2) { m_mur2=mur2; }
    inline void SetNorm(const double& norm) { m_norm=norm; }
    void SetCouplings(const MODEL::Coupling_Map& cpls);

    double AlphaQCD() const;
    double AlphaQED() const;

    inline int Mode() const    { return m_mode;    }
    inline int DRMode() const  { return m_drmode;  }
    inline int ColMode() const { return m_colmode; }

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

#define DECLARE_VIRTUALME2_GETTER(NAME,TAG)			   \
  DECLARE_GETTER(NAME,TAG,Virtual_ME2_Base,Process_Info);	   \
  void ATOOLS::Getter<Virtual_ME2_Base,Process_Info,NAME>::	   \
  PrintInfo(std::ostream &str,const size_t width) const		   \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
