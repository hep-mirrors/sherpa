#ifndef NLL_JetRate_H
#define NLL_JetRate_H

#include "Function_Base.H"
#include "MathTools.H"
#include "Gauss_Integrator.H"
#include "Fast_Function.H"
#include "NLL_Sudakov_Base.H"
#include "NLL_Branching_Probability_Base.H"
#include <map>

namespace SHERPA {

  typedef std::map<std::string,ATOOLS::Fast_Function *> Func_List;

  struct nll {
    enum code {
      Rate = 0,
      GqDg = 1,
      GgDg = 2,
      GfDf = 3,
      GstarDg  = 6,
      GqqgDg   = 7,
      IntGqDg  = 11,
      IntGgDg  = 12,
      IntGfDf  = 13,
      IntGqqgDg= 17,
      GqDgIntGgDg_GfDf= 21,
      GgDgIntGgDg     = 24,
      GgDgIntGfDf     = 25,
      GfDfIntGstarDg  = 26,
      GfDfIntGqqgDg   = 27,
      IntGgDgIntGgDg    = 34,
      IntGgDgIntGfDf    = 35,
      IntGfDfIntGstarDg = 36,
      IntGfDfIntGqqgDg  = 37,
      GqDgIntgg2_gggg_ggff_ffsg      = 46,
      GqDgIntffgg_gg2_ggff_gggg_ggtg = 47,

      GQDG     = 101,
      GgDG     = 102,
      GFDF     = 103,

      GqDG     = 104, //***

      IntGQDG     = 111,
      IntGgDG     = 112,
      IntGFDF     = 113,
      GQDGIntGgDG_GFDF = 121,
      GqDGIntGgDG_GfDf_GFDF =124, //***
      GQDGIntGgDG_GfDf_GFDF =125 //***
    };
  };

  class NLL_JetRate : public ATOOLS::Function_Base {
    nll::code m_mode;

    bool   m_table;
    int    m_massive;
    int    m_njet;
    double m_qmin,m_qmax;
    double m_qmin_b,m_qmax_b;

    static Func_List integrals;
    ATOOLS::Gauss_Integrator m_gauss;

    NLL_Sudakov_Base * p_deltaq, * p_deltag;
    NLL_Branching_Probability_Base * p_gammaq, * p_gammag, * p_gammaf;
    NLL_Sudakov_Base * p_deltaqm, * p_deltagm;
    NLL_Branching_Probability_Base * p_gammaqm, * p_gammafm;

    double Integrate(nll::code,double);
    double Integrate(nll::code,double,double);

    void InitIntegral(nll::code id, nll::code int_id, std::string name) ;
    void InitMassless();
    void InitMassive();
  public:
    NLL_JetRate(int njet ,double q0, double Q, int mass_flag);
    void SetSudakovs(NLL_Sudakov_Base *,NLL_Sudakov_Base *,NLL_Sudakov_Base * =0,NLL_Sudakov_Base * =0);
    void SetGammas(NLL_Branching_Probability_Base *,NLL_Branching_Probability_Base *,
		   NLL_Branching_Probability_Base *,NLL_Branching_Probability_Base * =0,
		   NLL_Branching_Probability_Base * =0);
    void Init();
    void Rates(double&,double&,double&,double&);
    double operator()(double);

    double SetQ0(double q0) { return m_qmin=q0; }
    double SetQ1(double q1) { return m_qmax=q1; }
    ~NLL_JetRate();
  };

}

#endif
