#ifndef NLL_Sudakov_H
#define NLL_Sudakov_H


#include "NLL_Sudakov_Base.H"
#include "NLL_Branching_Probability_Base.H"
#include "Message.H"
#include <vector>
#include <map>

namespace MODEL { class Running_AlphaS; }

namespace SHERPA {

  class NLL_JetRate;

  typedef std::map<ATOOLS::Flavour,NLL_Sudakov_Base*> Sudakov_Map;

  class NLL_Sudakov: public ATOOLS::Function_Base {
  private:

    int                     m_mode;
    MODEL::Running_AlphaS * p_runas;

    Sudakov_Map m_sudakovs;
    double      m_lambda, m_mu2, m_asmu, m_as_factor;

    void   PrepareMap();
    void   FixLambda2();    

    inline double Delta(const ATOOLS::Flavour &fl, double Q,double q)
    { return Delta(fl)(Q,q); }

  public:

    NLL_Sudakov(int mode, double q2max, double q2min, 
		MODEL::Running_AlphaS *as=0,double asfac=1.);

    ~NLL_Sudakov();

    NLL_Sudakov_Base & Delta(const ATOOLS::Flavour &);

  };

  /*!
    \file NLL_Sudakov.H
    \brief contains the class SHERPA::NLL_Sudakov
  */
  /*!
    \class NLL_Sudakov.H
    \brief this class provides the numerical values for the sudakov form factors 
           used in the merging procedure of parton shower and matrix elements

    The main purpose of this class is to provide sudakov form factors for
    the merging procedure of parton shower and matrix elements (cf. class 
    Cluster_Partons). Consequently the main routine is Delta(const ATOOLS::Flavour &)
    which returns the sudakov form factor for a given flavour.
  */
  /*!
    \fn NLL_Sudakov::NLL_Sudakov(double q2max, double q2min);
    \brief special constructor

    Constructor is called with a minimal and a maximal \f$q^2\f$. Within 
    the given domain sudakov form factors are initialized. Consequently a map 
    is filled. If neccessary lookup tables are precalculated.
  */
  /*!
    \fn NLL_Sudakov::~NLL_Sudakov();
    \brief destructor

    removes all sudakov objects from memory
  */
  /*!
    \fn void NLL_Sudakov::PrepareMap();
    \brief initializes a map with all (massless) sudakov form factors needed.

    For each Flavour (d-, u-, s-, c-, b-quark or anti-quark, and gluon) the corresponding
    sudakov form factor (cf. NLL_Single_Sudakov and NLL_Combined_Sudakov) are initialized
    and put into a map, allowing a fast access to a sudakov form factor object for any flavour.

    For purpose of completeness a NLL_Dummy_Sudakov (always one) is added to the map, 
    which will be returned for any flavour without a corresponding sudakov form 
    factor.

    In order to keep track of all sudakov objects inserted in the map, a list of
    unique NLL_Sudakov_Base objects is maintained, and will be used for a proper
    removal at the end.

    The default massless integrated splitting functions are obtained by 
    GammaQ_Lambda, GammaG_Lambda, and   GammaF_Lambda.
  */
  /*!
    \fn void NLL_Sudakov::PrepareMassiveMap();
    \brief initializes a map with all massive sudakov form factors.

    For each Flavour (d-, u-, s-, c-, b-quark or anti-quark, and gluon) the corresponding
    sudakov form factor (cf. NLL_Single_Sudakov and NLL_Combined_Sudakov) are initialized
    and put into a map, allowing a fast access to a sudakov form factor object for any 
    flavour. The difference to routine PrepareMap() comes from utilizing different integrated 
    splitting functions (cf. Krauss and Rodrigo hep-ph/0303038). Consequently different quark
    flavour can no longer be mapped to one common sudakov (as in the massless case) but
    have all to have individual sudakov objects instead.
    
    The massive integrated splitting functions are obtained by 
    GammaQ_Lambda_Massive, GammaG_Lambda_Massive, and   GammaF_Lambda_Massive.
  */
  /*!   
    \fn NLL_Sudakov_Base & NLL_Sudakov::Delta(const ATOOLS::Flavour &);
    \brief returns a NLL_Sudakov_Base object for a given flavour.

    This the main access method to NLL sudakov form factors. It can be 
    called after a map has been initialized via PrepareMap() or PrepareMassiveMap().
    After that each call of this routine returns the corresponding sudakov form
    factor in form of a NLL_Sudakov_Base object. 
    For not strong interacting particles usually an NLL_Dummy_Sudakov is returned, 
    giving a one for all possible scale combinations.
  */
  /*!
    \fn    double NLL_Sudakov::Delta(const ATOOLS::Flavour &, double ,double );
    \brief returns the result of a sudakov form factor for a given flavour and scales

    This is a wrapper around the routine Delta(const ATOOLS::Flavour &). It has basically 
    the same features, but is provided for convinience of the user.

    The following two lines are equivalent:
    \verbatim

      dg = sud.Delta(Flavour(kf::gluon))(Q,Q0);
      dg = sud.Delta(Flavour(kf::gluon),Q,Q0);
    \endverbatim
  */
  /*!
    \fn void   NLL_Sudakov::FixLambda2();    
    \brief sets up some input paramter for the calculation of the sudakov form factors

    This routine sets all parameter (correlated to the running coupling) used
    in sudakov form factors. In most cases the main parameter is \f$\Lambda_{\rm QCD}\f$.
    In some parameterization it may as well be a renomaization scale \f$\mu\f$ and a
    corresponding value for the strong coupling \f$\alpha_S(\mu)\f$.
  */
  /*!
    \fn void NLL_Sudakov::CheckSudakovs();
    \brief here some consitency checks are performed and some results for 
           selected paramters are returned.
  */
}
#endif
