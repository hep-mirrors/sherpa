#ifndef Spacelike_Sudakov_H
#define Spacelike_Sudakov_H

#include "Backward_Splitting_Group.H"
#include "Flavour.H"
#include "Knot.H"
#include "MathTools.H"

namespace PDF {
  class PDF_Base;
}


namespace APACIC {
  class Sudakov_Tools;

  // dicing z and t for spacelike (initial state) shower evolution
  class Spacelike_Sudakov : public Backward_Splitting_Group {
    // global parameters
    Sudakov_Tools * tools;
    PDF::PDF_Base * pdfa;
    double          pt2min, pt2max; 
    double          alphaSmax, alphaSfix; 
    double          lambda2, b, t0;
    double          xe, emin, s;

    double          alphaQED;
    double          pt2minQED;

    int             cpl_scheme;
    int             ordering_scheme;
    // mother
    double E2;
    double x,z,t,phi,pt2;
    double maxpt;
    APHYTOOLS::Flavour inflav;

    void Add(Splitting_Function * spl) ;
    double CrudeInt(double _zmin, double _zmax);
    
    void SelectOne() { selected->SelectOne();};
    void ProduceT();  
    bool Veto(Knot *); 
    bool CplVeto();
    bool PTVeto(Knot *);
    bool MassVeto();
    void UniformPhi() { phi =  2.*M_PI*AMATOOLS::ran.Get(); };
  public:
    Spacelike_Sudakov(PDF::PDF_Base *,Sudakov_Tools * _tools);
    ~Spacelike_Sudakov() {};
    bool Dice(Knot * mother,double sprime);

    static double Q02(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSMass()); 
    }; // !!!
    static double MinE2(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSMass()); 
    }; //!!!mass

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

  /*
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
#endif // Spacelicke_Sudakov_H









