#ifndef Timelike_Sudakov_H
#define Timelike_Sudakov_H

#include "Splitting_Group.H"
#include "Random.H"

namespace APACIC {
  class Sudakov_Tools;
  class Timelike_Kinematics;
  class Knot;

  class Timelike_Sudakov : public Splitting_Group {
    // global parameters
    Sudakov_Tools       * p_tools;
    Timelike_Kinematics * p_kin;
    double m_pt2min, m_pt2max; 
    double m_alphaSmax, m_alphaSfix; 
    double m_pt2minQED;
    double m_alphaQED;

    // switches:( get in constructor)
    int m_cpl_scheme; 
    int m_ordering_scheme; 
    int m_pt_scheme; 
    int m_mass_scheme;
    int m_width_scheme;
    int m_MEcorr_scheme;
    int m_zrange_scheme;
    int m_angle_scheme;
    int m_direct_photons;
  
    // mother
    double m_wa;
    double m_ta, m_tb, m_tc, m_t0;
    double m_z, m_phi, m_pt2;
    double m_maxpt;
    APHYTOOLS::Flavour m_inflav;

    // veto
    int m_last_veto;

    void CheckSplittings();
  public:
    Timelike_Sudakov(Timelike_Kinematics *);
    ~Timelike_Sudakov();

    bool Dice(Knot * mother,Knot * granny=0);
  private:
    double CrudeInt(double _zmin, double _zmax) ;

    void Add(Splitting_Function * spl);
    void SelectOne();

    void ProduceT();  

    bool Veto(Knot *); 

    bool CplVeto();
    bool AngleVeto(Knot *);
    bool MEVeto(Knot *);
    bool MassVeto();
    bool JetVeto(Knot *);

    void UniformPhi() { m_phi = 2.*M_PI*AMATOOLS::ran.Get(); }
  };


  /*!
    \file
    \brief contains the class APACIC::Timelike_Sudakov
  */

  /*!
    \class Timelike_Sudakov
    \brief Class for dicing z and t for timelike (final state) shower evolution
  */

  //-----------------------------------------------------------------------
  //-------------------- Constructors -------------------------------------
  //----------------------------------------------------------------------- 
  /*! 
    \fn Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics *);
    \brief Standard Constructor -

      loads of stuff to be removed ....
  */
  /*!
    \fn  Timelike_Sudakov::~Timelike_Sudakov();
    \brief Destructor
  */

  //-----------------------------------------------------------------------
  //-------------------- Dicing the next branch ---------------------------
  //----------------------------------------------------------------------- 
  /*! 
    \fn bool Timelike_Sudakov::Dice(Knot * mother,Knot * granny=0);
    \brief Tries to generate t and z for a mother ( returns "1" if successful),
       needs virtuality of granny and energy and flavour of mother
  */
  
  //-----------------------------------------------------------------------
  //-------------------- Handling pf splitting functions ------------------
  //----------------------------------------------------------------------- 
  /*! 
    \fn    void Timelike_Sudakov::Add(Splitting_Function * spl);
    \brief Initialises lists of splitting functions grouped by the flavour of the
           incoming leg A.
  */
  
  /*!
     \fn     double Timelike_Sudakov::CrudeInt(double _zmin, double _zmax) ;
     \brief Yields the crude integral of a group.
  */

     
  /*!
    \fn void Timelike_Sudakov::SelectOne();
    \brief Selects one specific mode for the splitting (to be called after CrudeInt)
  */

  //-----------------------------------------------------------------------
  //-------------------- Methods for dicing -------------------------------
  //----------------------------------------------------------------------- 
  /*! 
    \fn void Timelike_Sudakov::ProduceT();  
    \brief Starting with a t a new, a lower t will be produced and set.
    For quarks there's also the possibility to have q->q+photon. 
  */
  /*!
    \fn bool Timelike_Sudakov::Veto(Knot *);  
    \brief Organizes the vetoes to cure overshoot in approximate Sudakov formfactor
    used in ProduceT, returns a 0 or 1 if no veto or a veto is thrown.
    Calls the specific vetos after some trivial tests.
  */
  
  /*!
    \fn     bool  Timelike_Sudakov::CplVeto();
    \brief Corrects for overshooting in taking a crude alpha_S
  */

  /*!
    \fn     bool  Timelike_Sudakov::AngleVeto(Knot *);
    \brief Veto from imposed angular ordering on subsequent emissions
  */
  
  /*!
    \fn   bool Timelike_Sudakov::MEVeto(Knot *);
    \brief applies ME corrections to first correction

    Depending on ordering schemes,emissions get corrected on
    the ME kinematics of e+e- -> qqg
    if "first" branch perform ME - Correction for gluon radiation
    \verbatim
             (t',z') *      
               /     \      x_i = 2 E_i / sqrt(t')
            (t,z) *   \
             / \       \
            1   3       2
    \endverbatim
  */
  
  /*! 
    \fn bool Timelike_Sudakov::MassVeto();
    \brief checks for \e z being constient with physical range

    First veto from condition on cos(theta_bc) to be physical, condition:
    \f$ z_min < z < z_max\f$,   \f$z_min,max = 1/2 (1 -+ p/E)\f$
    Here, there's a constrained and an unconstrained choice.
    Second veto from modified splitting function for heavy quark a la
    Catani. 
  */
  
  /*!
    \fn bool Timelike_Sudakov::JetVeto(Knot *);  
    \brief considers only branching that do not constitute an independen jet
  */
    

  /*!
    \fn     void   Timelike_Sudakov::UniformPhi()
    \brief Uniform distribution of phi, angular correlations to be imposed.
  */


  /*!
    \fn void Timelike_Sudakov::CheckSplittings()
    \brief Perfroming consistency checks; for Internal use only
  */

} // namespace APACIC
#endif
