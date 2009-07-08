#ifndef APACIC_Showers_QED_Splitting_Functions_H
#define APACIC_Showers_QED_Splitting_Functions_H
//----------------------------------------------------------------------
//               QED_Splitting_Functions.H
//----------------------------------------------------------------------
#include "AddOns/Apacic++/Showers/Splitting_Function.H"

namespace APACIC {
  class Sudakov_Tools;

  // fermion to fermion + photon  splitting function
  class f_fp : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
    double          m_qsqr;
  public:
    f_fp(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour);
    f_fp(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double _zmin, double _zmax);
  };


  // fermion to photon + fermion splitting function (only used in Initial State Shower)
  class f_pf : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
    double          m_qsqr;
  public:
    f_pf(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour);
    f_pf(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double _zmin, double _zmax);
  };

  // photon to fermion + anti-fermion splitting function
  class p_ff : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha;
    double          m_qsqr;
  public:
    p_ff(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour);
    p_ff(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    double operator()(double z);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool masses);
    double CrudeInt(double _zmin, double _zmax);
  };


  /*! 
    \file
    \brief contains the classes APACIC::f_fp, APACIC::f_pf and APACIC::p_ff
  */

  /*!
    \class f_fp
    \brief represents fermion to fermion + photon  splitting function

    This class represents the splitting function for a \f$f \to f \; \gamma\f$
    splitting:
    \f[
       P_{f \to f \gamma} (z) = Q_f^2 \frac{1 + z^2}{1 -z}
    \f]
  */

  /*!
    \fn f_fp::f_fp(ATOOLS::Flavour fermionflavour);
    \brief Constructor, taking the fermion flavour as argument (obsolete)
  */

  /*!
    \fn f_fp::f_fp(ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    \brief Constructor, takng the fermion flavour and a class to
           deliver running alpha_qed 
  */

  /*!
    \class f_pf
    \brief represents fermion to photon + fermion splitting function 
         (only used in Initial State Shower)

    This class represents the splitting function for a \f$f \to \gamma \; f \f$
    splitting:
    \f[
       P_{f \to \gamma f} (z) = Q_f^2 \frac{1 + (1-z)^2}{z}
    \f]
  */
  /*!
    \fn f_pf::f_pf(ATOOLS::Flavour fermionflavour);
    \brief Constructor, taking the fermion flavour as argument (obsolete)
  */
  /*!
    \fn f_pf::f_pf(ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    \brief Constructor, takng the fermion flavour and a class to
           deliver running alpha_qed 
  */


  /*!
    \class p_ff
    \brief represents photon to fermion + anti-fermion splitting function

    This class represents the splitting function for a \f$\gamma \to f \; \bar f\f$
    splitting:
    \f[
       P_{\gamma \to f \bar f} (z) = Q_f^2 [ z^2 + (1 -z)^2 ]
    \f]
  */

  /*!
    \fn p_ff::p_ff(ATOOLS::Flavour fermionflavour);
    \brief Constructor, taking the fermion flavour as argument (obsolete)
  */
  /*!
    \fn p_ff::p_ff(ATOOLS::Flavour fermionflavour,Sudakov_Tools * _tools);
    \brief Constructor, takng the fermion flavour and a class to
           deliver running alpha_qed 
  */

} // namespace APACIC
#endif

