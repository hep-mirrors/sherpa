#ifndef timelike_sudakov_h
#define timelike_sudakov_h


#include "Sudakov_Tools.H"
#include "Splitting_Group.H"
#include "Vector.H"
#include "Flavour.H"
#include "Knot.H"
#include "MathTools.H"

namespace APACIC {
  class Timelike_Kinematics;

  // dicing z and t for timelike (final state) shower evolution
  class Timelike_Sudakov : public Splitting_Group {
    // global parameters
    Sudakov_Tools * tools;
    Timelike_Kinematics    * kin;
    double          pt2min, pt2max; 
    double          alphaSmax, alphaSfix; 
    double          pt2minQED;
    double          alphaQED;

    // switches:( get in constructor)
    int cpl_scheme; 
    int ordering_scheme; 
    int pt_scheme; 
    int mass_scheme;
    int width_scheme;
    int MEcorr_scheme;
    int zrange_scheme;
    int angle_scheme;
    int direct_photons;
  
    // mother
    double wa;
    double ta,tb,tc,t0;
    double z,phi,pt2;
    double maxpt;
    APHYTOOLS::Flavour inflav;

    void CheckSplittings();
  public:
    //-----------------------------------------------------------------------
    //-------------------- Constructors -------------------------------------
    //----------------------------------------------------------------------- 
    /*! 
      Standard Constructor -
      loads of stuff to be removed ....
    */
    Timelike_Sudakov(Timelike_Kinematics *);
    ~Timelike_Sudakov() {
      if (tools) delete tools;
    };

    //-----------------------------------------------------------------------
    //-------------------- Dicing the next branch ---------------------------
    //----------------------------------------------------------------------- 
    /*! 
       Tries to generate t and z for a mother ( returns "1" if successful),
       needs virtuality of granny and energy and flavour of mother
    */
    bool Dice(Knot * mother,Knot * granny=0);
  private:
    //-----------------------------------------------------------------------
    //-------------------- Handling pf splitting functions ------------------
    //----------------------------------------------------------------------- 
    /*! 
      Initialises lists of splitting functions grouped by the flavour of the
      incoming leg A.
    */
    void Add(Splitting_Function * spl) {
      for (SplFunIter iter(group);iter();++iter) {
	if (iter()->GetFlA()==spl->GetFlA()) {
	  iter()->Add(spl);
	  return ;
	}
      }
      group.Append(new Splitting_Group(spl));
      selected=spl;
    };
    /*!
      Yields the crude integral of a group.
    */
    double CrudeInt(double _zmin, double _zmax) {
      SplFunIter iter(group);
      for (;iter();++iter)
	if (iter()->GetFlA()==inflav) {selected=iter();break;}
      if (!iter()) {
	//cout<<"Timelike_Sudakov::CrudeInt : "<<endl;
	//cout<<"Warning : splitting function missing for "<<inflav<<endl;
	lastint = -1.;
      }
      return lastint=selected->CrudeInt(_zmin,_zmax);
    };        
    //! Selects one specific mode for the splitting (to be called after CrudeInt)
    void SelectOne() { selected->SelectOne(); };

    //-----------------------------------------------------------------------
    //-------------------- Methods for dicing -------------------------------
    //----------------------------------------------------------------------- 
    /*! 
      Starting with a t a new, a lower t will be produced and set.
      For quarks there's also the possibility to have q->q+photon. 
    */
    void ProduceT();  
    /*! 
      Organizes the vetoes to cure overshoot in approximate Sudakov formfactor
      used in ProduceT, returns a 0 or 1 if no veto or a veto is thrown.
      Calls the specific vetos after some trivial tests.
    */
    bool Veto(Knot *); 
    //! Corrects for overshooting in taking a crude alpha_S
    bool  CplVeto();
    //! Veto from imposed angular ordering on subsequent emissions
    bool  AngleVeto(Knot *);
    /*!
      Depending on ordering schemes,emissions get corrected on
      the ME kinematics of e+e- -> qqg
      if "first" branch perform ME - Correction for gluon radiation
             (t',z') *      
               /     \      x_i = 2 E_i / sqrt(t')
            (t,z) *   \
             / \       \
            1   3       2
    */
    bool MEVeto(Knot *);
    /*! 
      First veto from condition on cos(theta_bc) to be physical, condition:
      z_min < z < z_max,   z_min,max = 1/2 (1 -+ p/E)
      Here, there's a constrained and an unconstrained choice.
      Second veto from modified splitting function for heavy quark a la
      Catani. 
     */
    bool MassVeto();
    /*
      
    */
    bool JetVeto(Knot *);

    //! Uniform distribution of phi, angular correlations to be imposed.
    void   UniformPhi() { phi =  2.*M_PI*AMATOOLS::Ran.get(); };
  };
}
#endif





