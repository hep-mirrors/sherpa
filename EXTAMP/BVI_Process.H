#ifndef EXTAMP_BVI_Process_H
#define EXTAMP_BVI_Process_H

#include "EXTAMP/Process.H"

namespace PHASIC {
  class Virtual_ME2_Base;
  class Color_Correlated_ME2;
  class Process_Info;
  class KP_Terms;
}

namespace EXTAMP {

  class BVI_Process : public Process {

  public :

    BVI_Process(const PHASIC::Process_Info& pi,
		const double& vfrac,
		const int& subtraction_type);
    
    ~BVI_Process();

    double Partonic(const ATOOLS::Vec4D_Vector &p,
		    const int mode);

    void SetNLOMC(PDF::NLOMC_Base *const mc);

  private :

    /* Subtraction scheme: O: plain CS
                           1: Dire
			   2: modfied CS */
    int m_subtype;

    /* Evaluate virtuals only for a fraction m_vfrac of PS points */
    double m_vfrac;

    PHASIC::Virtual_ME2_Base*     p_loop_me;
    PHASIC::Color_Correlated_ME2* p_corr_me;
    PHASIC::KP_Terms*             p_kpterms;

    double Calc_V(const ATOOLS::Vec4D_Vector& p,
		  const double& born,
		  const double& mur) const;

    double Calc_I(const ATOOLS::Vec4D_Vector& p,
		  const double& mur) const;

    /* Can't be const because we need to store momentum fractions
       m_x0,m_x1,m_eta0,m_eta1 of KP_Terms for reweighting */
    double Calc_KP(const ATOOLS::Vec4D_Vector& p);
    double m_x0, m_x1, m_eta0, m_eta1;

    bool ApplyIDdipole() const;
    double Calc_HP(const ATOOLS::Vec4D_Vector& p);

    /* Used by PHASIC::Single_Process for reweighting op KP terms */
    double KPTerms(int mode, double scalefac2=1.0);


    /* Calc all terms explicitly dependent on renormalization scale
       for on-the-fly reweighting. First item: first derivative of all
       terms with respect to scale log, second item: second derivative
       of all terms with respect to scale log. */
    std::pair<double,double>
    Calc_ScaleDependenceTerms_I(const ATOOLS::Vec4D_Vector& p,
				const double& mur) const;
    std::pair<double,double>
    Calc_ScaleDependenceTerms_V(const ATOOLS::Vec4D_Vector& p,
				const double& B,
				const double& mur) const;

    /* epsilon^0, epsilon^{-1}, epsilon^{-2} coefficients of eq.
       (5.32) through (5.34) of arXiv:hep-ph/9605323, linearly
       combined according to (8.13) and (8.14) */
    static double Vi_eps0(const ATOOLS::Flavour& fl, int subtype);
    static double Vi_eps1(const ATOOLS::Flavour& fl);
    static double Vi_eps2(const ATOOLS::Flavour& fl);

    constexpr static double m_CF = 4./3.;
    constexpr static double m_CA = 3.0;
    constexpr static double m_TR = 1./2.;
    static double m_NF;

    /* beta0 = 11/3 CA - 4/3 TR nf*/
    double m_beta0;
    
    static double Ti2(const ATOOLS::Flavour& fl);

  };

}

#endif
