#ifndef Spacelike_Sudakov_H
#define Spacelike_Sudakov_H

#include "Backward_Splitting_Group.H"
#include "Random.H"

namespace PDF {
  class PDF_Base;
}

namespace APACIC {
  class Sudakov_Tools;
  class Spacelike_Kinematics;
  class Knot;

  // dicing z and t for spacelike (initial state) shower evolution
  class Spacelike_Sudakov : public Backward_Splitting_Group {
    // global parameters
    Sudakov_Tools * p_tools;
    Spacelike_Kinematics * p_kin;
    PDF::PDF_Base * p_pdfa;
    double          m_pt2min, m_pt2max; 
    double          m_alphaSmax, m_alphaSfix; 
    double          m_lambda2, m_b, m_t0;
    double          m_xe, m_emin, m_s;

    double          m_alphaQED;
    double          m_pt2minQED;

    int             m_cpl_scheme;
    int             m_ordering_scheme;
    int             m_jetveto_scheme;
    int             m_last_veto;
    double          m_qjet;
    double          m_pdf_fac;
    // mother
    double m_E2;
    double m_x, m_z, m_t, m_phi,m_pt2;
    double m_maxpt;
    APHYTOOLS::Flavour m_inflav;

    void Add(Splitting_Function * spl) ;
    double CrudeInt(double _zmin, double _zmax);
    
    void SelectOne() { p_selected->SelectOne();};
    void ProduceT();  
    bool Veto(Knot *,bool,int & extra_pdf); 
    bool CplVeto();
    bool PTVeto(Knot *);
    bool JetVeto(Knot *);
    bool MassVeto(int);
    void UniformPhi() { m_phi =  2.*M_PI*AMATOOLS::ran.Get(); };
  public:
    Spacelike_Sudakov(PDF::PDF_Base *,Sudakov_Tools *, Spacelike_Kinematics *,
		      double,AORGTOOLS::Data_Read *); 
    ~Spacelike_Sudakov() {};
    bool Dice(Knot * mother,double sprime,bool jetveto, int & extra_pdf);

    static double Q02(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSMass()); 
    } // !!!
    static double MinE2(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSMass()); 
    } //!!!mass

  };

  /*!
    \file
    \brief contains the class APACIC::Spacelike_Sudakov
  */

  /*!
    \class Spacelike_Sudakov
    \brief provides all necessary functions to calculate sudakovs, and 
    perform corresponding shower evolution.
  */

  /*!
    \fn    void Spacelike_Sudakov::Add(Splitting_Function * spl) 
    \brief Initialises lists of splitting functions grouped by the flavour of leg B.
           Note that it is a backward evolution!
  */

  /*!
    \fn double Spacelike_Sudakov::CrudeInt(double _zmin, double _zmax) {
    \brief Yields the crude integral of a group.
  */
    
  /*!
    \fn void Spacelike_Sudakov::SelectOne() { selected->SelectOne();};
    \brief Selects on particluar splitting mode out of a group.
  */
  
  /*!
    \fn void Spacelike_Sudakov::ProduceT();  
    \brief  Starting with a t a new, "lower" virtuality will be produced and set.
  */
    
  /*! 
    \fn bool Spacelike_Sudakov::Veto(Knot *); 
    \brief Organizes the vetoes to cure overshoot in approximate Sudakov formfactor
      used in ProduceT, returns a 0 or 1 if no veto or a veto is thrown.
      Calls the specific vetos after some trivial tests.
  */

  /*!
    \fn bool  Spacelike_Sudakov::CplVeto();
    \brief     Corrects for overshot in taking a crude coupling
  */

  /*!
    \fn    bool  Spacelike_Sudakov::PTVeto(Knot *);
    \brief Veto from imposed pt ordering of previous branches. Note : pt = (1-z) t
  */

  /*! 
    \fn    bool Spacelike_Sudakov::MassVeto();
    \brief Check z ranges, weight from pdfs as well. 
  */

  /*!
    \fn    void   Spacelike_Sudakov::UniformPhi();
    \brief Returns an azimuthal angle uniformly distributed.
  */

  /*! 
    \fn  Spacelike_Sudakov::Spacelike_Sudakov(PDF::PDF_Base *,Sudakov_Tools * _tools);
    \brief Standard Constructor -
      loads of stuff to be removed ....
  */
  
  /*!
    \fn Spacelike_Sudakov::~Spacelike_Sudakov();
    \brief Destructor
  */

  /*! 
    \fn bool Spacelike_Sudakov::Dice(Knot * mother,double sprime);
    \brief  Tries to generate t and z  ( returns "1" if successful)
      needs virtuality of daughter (remember : Backwards evolution !)
      and energy and flavour of mother and the actual sprime.
  */
  
  /*!
    \fn static double Spacelike_Sudakov::Q02(APHYTOOLS::Flavour const & flav);
    \brief returns flavour dependent Q0
  */

  /*!
    \fn static double Spacelike_Sudakov::MinE2(APHYTOOLS::Flavour const & flav)
    \brief returns flavour dependent Emin
  */



} // namespace APACIC
#endif // Spacelike_Sudakov_H
