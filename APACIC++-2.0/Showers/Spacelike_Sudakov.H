#ifndef Spacelike_Sudakov_h
#define Spacelike_Sudakov_h


#include "Timelike_Sudakov.H"
#include "Sudakov_Tools.H"
#include "Backward_Splitting_Group.H"
#include "Splitting_Group.H"
#include "Vector.H"
#include "Flavour.H"
#include "Knot.H"
#include "MathTools.H"


namespace APACIC {
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

    /*! 
      Initialises lists of splitting functions grouped by the flavour of leg B.
      Note that it is a backward evolution!
    */
    void Add(Splitting_Function * spl) {
      for (SplFunIter iter(group);iter();++iter) {
	if (iter()->GetFlB()==spl->GetFlB()) {
	  iter()->Add(spl);
	  return ;
	}
      }
      group.Append(new Backward_Splitting_Group(spl,pdf));
      selected=spl;
    };
    /*!
      Yields the crude integral of a group.
    */
    double CrudeInt(double _zmin, double _zmax) {
      SplFunIter iter(group);
      for (;iter();++iter)
	if (iter()->GetFlB()==inflav) { selected=iter(); break; }
      if (!iter()) std::cout<<"warning "<<inflav<<" not implemented in Splitting_Group"<<std::endl;
      return lastint = selected->CrudeInt(_zmin,_zmax);
    };    
    
    /*!
      Selects on particluar splitting mode out of a group.
    */
    void SelectOne() { selected->SelectOne();};
    //! Starting with a t a new, "lower" virtuality will be produced and set.
    void ProduceT();  
    /*! 
      Organizes the vetoes to cure overshoot in approximate Sudakov formfactor
      used in ProduceT, returns a 0 or 1 if no veto or a veto is thrown.
      Calls the specific vetos after some trivial tests.
    */
    bool Veto(Knot *); 
    //! Corrects for overshot in taking a crude coupling
    bool  CplVeto();
    //! Veto from imposed pt ordering of previous branches. Note : pt = (1-z) t
    bool  PTVeto(Knot *);
    /*! 
      Check z ranges, weight from pdfs as well. 
     */
    bool MassVeto();
    /*!
      Returns an azimuthal angle uniformly distributed.
    */
    void   UniformPhi() { phi =  2.*M_PI*AMATOOLS::Ran.get(); };
  public:
    /*! 
      Standard Constructor -
      loads of stuff to be removed ....
    */
    Spacelike_Sudakov(PDF::PDF_Base *,Sudakov_Tools * _tools);
    ~Spacelike_Sudakov() {};
    /*! 
      Tries to generate t and z  ( returns "1" if successful)
      needs virtuality of daughter (remember : Backwards evolution !)
      and energy and flavour of mother and the actual sprime.
    */
    bool Dice(Knot * mother,double sprime);

    static double Q02(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSmass()); 
    }; // !!!
    static double MinE2(APHYTOOLS::Flavour const & flav) { 
      return AMATOOLS::sqr(flav.PSmass()); 
    }; //!!!mass

  };

} // namespace APACIC
#endif // Spacelicke_Sudakov_H









