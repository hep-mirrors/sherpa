#ifndef Final_State_Shower_H
#define Final_State_Shower_H

#include "Timelike_Sudakov.H"
#include "Timelike_Kinematics.H"
#include "Tree.H"

#include "Blob_List.H"
#include "Blob.H"
#include "Parton_List.H"
#include "Run_Parameter.H"

namespace APACIC {
  /*!
    Does the final state parton shower evolution of a given particle ensemble.
    It will be called from outside via 
    1. PerformShowers
    2. FirstTimelikeFromSpacelike

    I think it should contain also all possible F.S. particle decays, but this 
    is something I still have to work on.
  */
  class Final_State_Shower {
    Timelike_Sudakov       * sud;
    Timelike_Kinematics    * kin;
    std::vector<Knot*> ini_partons;
  public:
    //-----------------------------------------------------------------------
    //--------------------------- Constructors ------------------------------
    //----------------------------------------------------------------------- 
    //! standard constructor
    Final_State_Shower();
    ~Final_State_Shower();

    //-----------------------------------------------------------------------
    //----------------------- Performing the Shower -------------------------
    //----------------------------------------------------------------------- 
    /*!
      Performs the FS shower on a given tree, i.e. starting from its root
      via InitializeJets(Tree,Tree->GetRoot) and does the kinematics afterwards.
      Return values are 1 complete success, 0 Kinematics didn't work out, 
      3 looks (almost) fine, but njet_ini is reduced (due to virtualities)!
    */
    int PerformShower(Tree *,bool);
    /*!
      Tries to initialize new jet system emerging from a timelike particle emitted from 
      a spacelike shower. If the Timelike_Sudakov yields a suitable virtuality t>t0, 
      daughters are initialized and the jet is evolved with EvolveJet. If the jet evolution
      was not successful a new t is generated. 
      NOTE : So far, only some vetos (coherence etc) have been implemented !!!
      Otherwise, the particle is set on its "mass shell" according to the definitions of 
      the Final_State_Shower and Timelike_Sudakov.
    */
    void FirstTimelikeFromSpacelike(Tree *,Knot *,bool);
    //-----------------------------------------------------------------------
    //---------------------------- After the Shower -------------------------
    //----------------------------------------------------------------------- 
    /*!
      Sets the colurs for the offsprings.
     */
    bool SetColours(Knot *);
    /*!
      After the shower is complete the outgoing knots of the tree are put into
      the parton list. Moreover, blobs are filled.
    */
    void ExtractPartons(Knot *,APHYTOOLS::Blob *,
			APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    // a simpler version (collects only final state partons)
    void ExtractPartons(Knot *,APHYTOOLS::Parton_List *);

    //-----------------------------------------------------------------------
    //---------------------------- Helpers ----------------------------------
    //----------------------------------------------------------------------- 
    /*!
      Test the Final State Shower with two toy-jets.
      Process will be   gamma* -> q bar q    at a fixed
      energy scale. 
    */
    bool TestShower(Tree *) ;
    //! check momentum conservation (adds daughter momenta recursively)
    AMATOOLS::vec4d GetMomentum(Knot * mo, int & number);
    /*!
      Outputs the tree.
    */
    void OutputTree(Tree *);

  private:
    //-----------------------------------------------------------------------
    //------------------- Initialisation of the Shower ----------------------
    //----------------------------------------------------------------------- 
    /*!
      Tries to initialize the jet system emerging from an incoming knot
      recursively along the following steps: 
      First, if any of the two daughters can not decay, i.e. is an internal line of
      a ME, InitializeJets is called with this daughter as an argument.
      If any of the two daughters can decay, the mothers branch is to be filled and
      the parton system produced by the decaying daughter(s) has to be evolved.
      This is repeated until either a fit has been found, then InitializeJets was 
      successful and a 1 is returned. If the evolution of the daughters was not successful,
      the initial kinematics get restored and the daughters are diced down. 
      In case no kinematically allowed daughter system can be reconstructed,
      a 0 is returned.      
      New feature in case the now massive daughters are not anymore separated to form 
      individual jet, a 3 is returned. (need new ME Kinematics)
     */
    int InitializeJets(Tree *,Knot *);

    /*!
      Distributes the masses of daughters according to a Breit-Wigner distribution
      in case their flavours are treated as unstable. The idea is to use the tout's
      of the daughters to reshuffle the momenta. In case a smearing occured
      a 1 is returned.
    */
    bool SmearDaughters(Knot * mo);

    /*! 
       Initialise two toy-jets - mainly used for self-tests of
       the Final_State_Shower etc..
    */
    void InitTwojetTree(Tree *,double);

    //-----------------------------------------------------------------------
    //------------------------ Evolution of the Shower ----------------------
    //----------------------------------------------------------------------- 
    /*!
      Evolves the parton system radiated off an given knot and returns
      1 when successful (0 if not). 

      This evolution is done recursively, where each step proceeds as follows :

      If the given mother knot is inactive, its daughters are deleted and the 
      evolution is succesful for this parton (return 1).
      Otherwise, if the branch gets filled, i.e. if daughters for the active
      mother with given t and z can exist after all, they get evolved further.
      If this works out a 1 is returned, otherwise the daughters and the mother
      z (which might have been altered during the evolution) are restored for
      another try, which dices down the t of the daughters. 
      This is repeated until suitable daughters are found or both daughters are inactive 
      with no fitting kinematics constructed. In the latter case, the evolution is not 
      successful (return value 0).
    */
    bool EvolveJet(Tree *,Knot *);

    /*!
      Try to branch a mother knot with given t and z :
      One of the daughter gets diced down more and the daughters system is
      checked for a kinematic fit. Then 1 is returned when a suitable branch 
      has been found, i.e. kinematically allowed daughters. Note, 
      that depending on the integer flag, either only the mothers z (i.e. the energy 
      fractions of the daughters) or the daughters momenta (if the daughters are from a 
      ME) get shuffled and eventually have to be restored.
    */
    bool FillBranch(Tree *,Knot *,int);

    //-----------------------------------------------------------------------
    //------------------- Service for the Branchings ------------------------
    //----------------------------------------------------------------------- 
    /*!
      New daughter to dice {0 or 1} for {d1 or d2}
      is chosen according to who's farer off-shell.
      Constraints on daughters apply if one of them
      is already inactive, i.e. its t is set on t_out.
     */
    bool ChooseDaughter(Knot *);

    /*!
      Initializes new daughters. If knots do not exist so far,
      they are initialized including there mother-daughter relation
      and the flavour they carry. If they are already initialised,
      their status is reset to "active", and their t and E^2 are reset
    */
    void InitDaughters(Tree *,Knot *,APHYTOOLS::Flavour *,bool);

    /*!
      Calls to reset the left and the right knot of the mother
    */
    void ResetDaughters(Knot *);

    /*!
      Deletes recursively the daughters and next generations of 
      each knot, partons are set inactive (leading to empty entries
      in the partonlist), they status of the knots is set to inactive,
      their t to zero and their daughters get deleted.
     */
    void Reset(Knot *);

    // make shure number of jet from hard event is not reduced!!!
    bool ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2);
    // make shure number of jet from hard event is not reduced using ALL ini_partons
    bool ExtraJetCheck();
  };

}
#endif



