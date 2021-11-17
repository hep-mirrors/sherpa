#ifndef HADRONS_Current_Library_VA_0_PP_H
#define HADRONS_Current_Library_VA_0_PP_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {    
  class VA_0_PP : public Current_Base {
    class FormFactor {
    protected:
      // a list of resonances and the prefactors that we use when combinung them
      // in the overall form factor calculation.
      std::list<std::pair<ResonanceFlavour, double> > m_Vresonances, m_Sresonances;
      double  m_Vresonance_norm, m_Sresonance_norm;
      virtual void FillResonances(const GeneralModel & model)=0;
    public:
      FormFactor(const GeneralModel & model) :
	m_Vresonance_norm(0.), m_Sresonance_norm(0.) {}
      virtual Complex F_V(const double & q2)=0;
      virtual Complex F_S(const double & q2)=0;
    };

  private:
    FormFactor * p_ff; // ff model
    double       m_global;   // V_CKM * Clebsch Gordon
    double       m_DeltaM2;  // mass difference of the two mesons
  public:
    VA_0_PP(const ME_Parameters &parameters,const std::string & name);
    ~VA_0_PP() { if (p_ff) { delete p_ff; p_ff = 0; } }
    virtual void SetModelParameters( struct GeneralModel _md ) {}
    virtual void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);

    /////////////////////////////////////////////////////////////////
    // Kuehn-Santamaria form factors
    /////////////////////////////////////////////////////////////////
    class FF_KS_pipi : public FormFactor {
      void FillResonances(const GeneralModel & model);
    public:
      FF_KS_pipi(const GeneralModel & model);
      Complex F_V(const double & q2);
      Complex F_S(const double & q2) { return Complex(0.,0.); }
    };

    class FF_KS_Kpi : public FormFactor {
      void FillResonances(const GeneralModel & model);
    public:
      FF_KS_Kpi(const GeneralModel & model);
      Complex F_V(const double & q2);
      Complex F_S(const double & q2);
    };

    /////////////////////////////////////////////////////////////////
    // Resonance Chiral Theory form factors
    /////////////////////////////////////////////////////////////////
    class FF_RChT_pipi : public FormFactor {
      double  m_fpi, m_m2_pi, m_m2_K;
      void    FillResonances(const GeneralModel & model);
      Complex A_Loop( double x, double y );
    public:
      FF_RChT_pipi(const GeneralModel & model);
      Complex F_V(const double & q2);
      Complex F_S(const double & q2) { return Complex(0.,0.); }
    };    

    class FF_RChT_Kpi : public FormFactor {
      double  m_m2_pi, m_m2_K, m_m2_eta, m_DeltaM2, m_SigmaM2, m_Sigma_KP, m_Delta_KP;
      double  m_piK_threshold, m_etaK_threshold;
      double  m_fpi, m_mK2, m_fpi2, m_mu2, m_mPi2, m_mEta2, m_MK02, m_cd, m_cm, q2;
      void    FillResonances(const GeneralModel & model);
    public:
      FF_RChT_Kpi(const GeneralModel & model);
      Complex F_V(const double & q2, double s);
      Complex F_S(const double & q2);
    };  
    /////////////////////////////////////////////////////////////////
    // Resonance Chiral Theory form factors
    /////////////////////////////////////////////////////////////////
    class FF_RChT2_pipi : public FormFactor {
      double  m_fpi, m_m2_pi, m_m2_K;
      void    FillResonances(const GeneralModel & model);
      Complex A_Loop( double x, double y );
    public:
      FF_RChT2_pipi(const GeneralModel & model);
      Complex F_V(const double & q2);
      Complex F_S(const double & q2) { return Complex(0.,0.);
      Complex CalcA(double m_m);}
    };    

    class FF_RChT2_Kpi : public FormFactor {
      double  m_m2_pi, m_m2_K, m_m2_eta, m_DeltaM2, m_SigmaM2;
      double  m_piK_threshold, m_etaK_threshold;
      double  m_fpi, m_fpi2, m_mu2, m_cd, m_cm;
      void    FillResonances(const GeneralModel & model);
    public:
      FF_RChT2_Kpi(const GeneralModel & model);
      Complex F_V(const double & q2, double s);
      Complex F_S(const double & q2);
      Complex CalcA(double m_m);};

    
    /////////////////////////////////////////////////////////////////
    // Resonance Gounaris-Sakurai Chiral Theory form factors
    /////////////////////////////////////////////////////////////////
    class FF_GS_pipi : public FormFactor {
      double  m_fpi, m_m2_pi, m_m2_K, q2;
      void    FillResonances(const GeneralModel & model);
      Complex A_Loop( double x, double y );
    public:
      Complex FF_RChT2_pipi(const GeneralModel & model);
      Complex F_V(const double & q2);
      Complex F_S(const double & q2) { return Complex(0.,0.); }
    };          
  };
}
#endif
