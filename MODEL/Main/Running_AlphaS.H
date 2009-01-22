#ifndef MODEL_Main_Running_AlphaS_H
#define MODEL_Main_Running_AlphaS_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"

namespace MODEL {
  class AS_Data {
  public:
    int    nF;
    double beta0,beta1,beta2,beta3;
    double mass2_low,mass2_up,as_low,as_up;
    double lambda2;
  };

  class AS_Sorting_Criterion {
  public :
    bool operator()(const AS_Data * data1,const AS_Data * data2) 
    {
      return (data1->mass2_up<data2->mass2_up);
    }
  };

  typedef std::set<AS_Data *,AS_Sorting_Criterion> AS_Data_set; 
  typedef AS_Data_set::iterator                    asdata_iterator; 

  class Running_AlphaS : public ATOOLS::Function_Base {
  protected:
    static double s_Nc, s_CA, s_CF, s_TF, s_zeta2, s_zeta3;

    double        m_asMZ, m_M2Z;
    int           m_order;
    double        m_fac, m_ascut, m_cutq2;

    AS_Data_set     m_asdata;
    asdata_iterator m_asditer;

    void FixASData();
    void FixLambda(AS_Data * region,
		   const double scale,const double asvalue) const;
    void FixSoftEnd();

    AS_Data *    FindRegion(const double q2) const;
    const double Beta0(const int nF) const;
    const double Beta1(const int nF) const;
    const double Beta2(const int nF) const;
    const double Beta3(const int nF) const;

    const double LambdaForm(const double q2,const AS_Data * asregion,
			    const int order) const;
    const double MatchingFactorDown(const AS_Data * up,const AS_Data * down,
				    const double refmass=-1.) const;
    const double MatchingFactorUp(const AS_Data * down,const AS_Data * up,
				  const double refmass=-1.) const;
  public:
    Running_AlphaS(const double as_MZ,const double m2_MZ,
		   const int order=0, const double fac=1.);
    ~Running_AlphaS();

    const double operator()(const double q2) const; 
    const double Inverse(const double as) const;
    const double AlphaS(const double q2) const;


    const int    Order() const { return m_order; }
    const double Beta0(const double q2) const;
    const double Beta1(const double q2) const;
    const double Beta2(const double q2) const;
    const double Beta3(const double q2) const;
    const int    Nf(const double q2) const;

    const double CutQ2() const;

    void   SelfTest() const;
  };



  extern Running_AlphaS * as;

  //! QCD SM Running of AlphaS
  /*!
    \class Running_AlphaS
    \brief The class for the (running) strong coupling constant.

    This is an implementation of the
    <A HREF="http://www-spires.dur.ac.uk/cgi-bin/spiface/find/hep/www?eprint=hep-ph/9706430">
    strong coupling constant by K.~G.~Chetyrkin, B.~A.~Kniehl and M.~Steinhauser</A> with the
    <A HREF="http://www-spires.dur.ac.uk/cgi-bin/spiface/hep/www?eprint=hep-ph/9701390">
    coefficients of the \f$\beta\f$-function by T.~van Ritbergen, J.~A.~M.~Vermaseren and 
    S.~A.~Larin </A>

    The following points should be noted:
    - In the two papers two different conventions are used for the coefficients of the 
      \f$beta\f$-function, which result in a relative factor of 4.  We use the
      definition of T.~van Ritbergen, J.~A.~M.~Vermaseren and S.~A.~Larin and adjust
      this in the code (with a comment).
    - At the moment the running at up to three loops is tested, the four loop
      running remains to be checked. 
    - Thresholds are fixed to the onshell masses of the corresponding quark.
    - In order to avoid the Landau Pole, the running is stopped as soon as 
      \f$alpha_s\f$ equals \f$4\pi\f$. (Below that scale \f$alpha_s\f$ is set 
      to constant.)
  */
  /*!
    \var int Running_AlphaS::m_order
    The order at which running \f$alpha_s\f$ is evaluated, e.g. m_order = 0 refers to one-loop,
    m_order = 3 to four-loop running.
  */



  /*!
    \var int Running_AlphaS::m_mzset
    The place of the Z-pole in the list of thresholds.
  */
  /*!
    \var double Running_AlphaS::m_m2_MZ
    The mass squared of the Z-boson.
  */
  /*!
    \var double Running_AlphaS::m_as_MZ
    \f$\alpha_s\f$ at the Z-pole.
  */
  /*!
    \var double Running_AlphaS::Running_AlphaS::m_CF
    The Casimir operator of QCD in the fundamental representation, \f$C_F=4/3\f$.
  */
  /*!
    \var double Running_AlphaS::Running_AlphaS::m_CA
    The Casimir operator of QCD in the adjoint representation, \f$C_A=3\f$.
  */
  /*!
    \var AsDataSet * Running_AlphaS::p_thresh
    Values of $\alpha_s$ and the \f[\beta_i\f] at the various thresholds.
  */
  /*!
    \fn double Running_AlphaS::Beta0(const double)
    Returns \f[\beta_0\f] at a given scale.
  */
  /*!
    \fn double Running_AlphaS::Beta0(const int)
    Returns \f[\beta_0 = (33 - 2n_f)/12\f] for a given number of active flavours.
  */
  /*!
    \fn double Running_AlphaS::Beta1(const int)
    Returns \f[\beta_1 = (306 - 38n_f)/48\f] for a given number of active flavours.
  */
  /*!
    \fn double Running_AlphaS::Beta2(const int)
    Returns \f[\beta_2 = (2857/2 - 5033/18n_f + 325/54n_f^2)/64\f] for a given number of active flavours.
  */
  /*!
    \fn double Running_AlphaS::Beta3(const int)
    Returns \f[\beta_3 = \left[ (149753/6 + 3564\zeta(3)) -
		                (1078361/162 + 6508/27\zeta(3))n_f +
		                (50065/162 + 6472/81\zeta(3))n_f^2 +
		                (1093/729)n_f^3\right]/256\f] for a given number of active flavours.
  */
  /*!
    \fn double Running_AlphaS::Lambda2(const int)
    Fills in the b-coefficients, \f[b_{1,2,3} = \beta_{1,2,3}/\beta_0\,,\f]
    and caluclates a suitable \f$\Lambda_{\rm QCD}^2\f$ which is then returned.
  */
  /*!
    \fn double Running_AlphaS::ZetaOS2(const double,const double,const double,const int)
    \todo Documentation for this method.
  */
  /*!
    \fn double Running_AlphaS::InvZetaOS2(const double,const double,const double,const int)
    \todo Documentation for this method.
  */
  /*!
    \fn double Running_AlphaS::AlphaSLam(const double,const int)
    \todo Documentation for this method.
  */
  /*!
    \fn void Running_AlphaS::ContinueAlphaS(int &)
    \todo Documentation for this method.
  */
  /*!
    \fn Running_AlphaS::Running_AlphaS(double,double)
    The constructor of Running_AlphaS fills in and sorts the thresholds according to their mass.
    \todo Documentation for this method.
  */
  /*!
    \fn Running_AlphaS::~Running_AlphaS()
    \todo Documentation for this method.
  */
  /*!
    \fn double Running_AlphaS::operator()(double) 
    Returns the runnnig \f$alpha_s\f$ at a given scale.
  */
  /*!
    \fn double Running_AlphaS::AlphaS(const double)
    Returns the value of the Running_AlphaS::operator()(double).
  */
  /*!
    \fn int Running_AlphaS::Nf(const double)
    Returns the number of active flavours at a given scale.
  */
  /*!
    \fn void Running_AlphaS::SelfTest()
    Prints out the results of a consistency check.
  */
}



#endif
