#ifndef Final_State_Shower_H
#define Final_State_Shower_H

#include <vector>

#include "Data_Read.H"
#include "Vector.H"
#include "Blob.H"
#include "Blob_List.H"
#include "Parton_List.H"

namespace APHYTOOLS {
  class Flavour;
}

namespace APACIC {
  class Timelike_Sudakov;
  class Timelike_Kinematics;
  class Knot;
  class Tree;

  class Final_State_Shower {
    double                   m_pt2min;
    Timelike_Sudakov       * p_sud;
    Timelike_Kinematics    * p_kin;
    std::vector<Knot*> m_ini_partons;
  public:
    Final_State_Shower(AORGTOOLS::Data_Read *);
    ~Final_State_Shower();

    int PerformShower(Tree *,bool);
    void FirstTimelikeFromSpacelike(Tree *,Knot *,bool,double,double);
    bool SetAllColours(Knot *);
    static bool SetColours(Knot *,Timelike_Kinematics *);
    static void EstablishRelations(Knot * mo, Knot * d1,Knot * d2);


    void ExtractPartons(Knot *,APHYTOOLS::Blob *,
			APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List * = NULL);
    void ExtractPartons(Knot *,APHYTOOLS::Parton_List * = NULL);

    bool TestShower(Tree *) ;
    AMATOOLS::Vec4D GetMomentum(Knot * mo, int & number);
    void OutputTree(Tree *);
    double PT2Min() { return m_pt2min; }
  private:
    int InitializeJets(Tree *,Knot *,int =0);
    bool SmearDaughters(Knot * mo);

    void InitTwojetTree(Tree *,double);
    bool EvolveJet(Tree *,Knot *);
    bool FillBranch(Tree *,Knot *,int);

    bool ChooseDaughter(Knot *);

    void InitDaughters(Tree *,Knot *,APHYTOOLS::Flavour *,bool);
    void ResetDaughters(Knot *);

    void Reset(Knot *);

    bool ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2);
    bool ExtraJetCheck();
  };


  /*!
    \file
    \brief contains the class APACIC::Final_State_Shower
  */

  /*!
    \class Final_State_Shower
    \brief    Does the final state parton shower evolution of a given particle ensemble.

    Does the final state parton shower evolution of a given particle ensemble.
    It will be called from outside via 
     -# PerformShower() or
     -# FirstTimelikeFromSpacelike()

    I think it should contain also all possible F.S. particle decays, but this 
    is something I still have to work on.
  */

  //-----------------------------------------------------------------------
  //--------------------------- Constructors ------------------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn Final_State_Shower::Final_State_Shower();
    \brief standard constructor
  */

  /*!
    \fn Final_State_Shower::~Final_State_Shower();
    \brief destructor
  */

  //-----------------------------------------------------------------------
  //----------------------- Performing the Shower -------------------------
  //----------------------------------------------------------------------- 
  /*!
     \fn int Final_State_Shower::PerformShower(Tree *,bool);
     \brief Performs the FS shower on a given tree, i.e. starting from its root
      via InitializeJets(Tree,Tree->GetRoot) and does the kinematics afterwards.
      Return values are 1 complete success, 0 Kinematics didn't work out, 
      3 looks (almost) fine, but njet_ini is reduced (due to virtualities)!
  */
  
  /*!
    \fn    void Final_State_Shower::FirstTimelikeFromSpacelike(Tree *,Knot *,bool,double,double);
    \brief Tries to initialize new jet system emerging from a timelike particle emitted from 
      a spacelike shower. If the Timelike_Sudakov yields a suitable virtuality t>t0, 
      daughters are initialized and the jet is evolved with EvolveJet. If the jet evolution
      was not successful a new t is generated. 
      Otherwise, the particle is set on its "mass shell" according to the definitions of 
      the Final_State_Shower and Timelike_Sudakov.
  */
  //-----------------------------------------------------------------------
  //---------------------------- After the Shower -------------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn bool Final_State_Shower::SetColours(Knot *);
    \brief Sets the colurs for the offsprings.
  */
  
  /*
    \fn void Final_State_Shower::ExtractPartons(Knot *,APHYTOOLS::Parton_List *);
    \brief a simpler version of the method below (collects only final state partons)
  */

  /*!
    \fn void Final_State_Shower::ExtractPartons(Knot *,APHYTOOLS::Blob *,
		APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *);
    \brief After the shower is complete the outgoing knots of the tree are put into
    the parton list. Moreover, blobs are filled.
  */

  //-----------------------------------------------------------------------
  //---------------------------- Helpers ----------------------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn bool Final_State_Shower::TestShower(Tree *);
    \brief Test the Final State Shower with two toy-jets.
    Process will be   gamma* -> q bar q    at a fixed
    energy scale. 
  */

  /*!
    \fn AMATOOLS::Vec4D Final_State_Shower::GetMomentum(Knot * mo, int & number);
    \brief check momentum conservation (adds daughter momenta recursively)
  */

  /*!
    \fn void Final_State_Shower::OutputTree(Tree *);
    \brief Outputs the tree.
  */
  
  //-----------------------------------------------------------------------
  //------------------- Initialisation of the Shower ----------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn int Final_State_Shower::InitializeJets(Tree *,Knot *);
    \brief Tries to initialize the jet system

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
    
  /*!
    \fn bool Final_State_Shower::SmearDaughters(Knot * mo);
    \brief   Distributes the masses of daughters according to a Breit-Wigner distribution
    in case their flavours are treated as unstable. The idea is to use the tout's
    of the daughters to reshuffle the momenta. In case a smearing occured
    a 1 is returned.
  */
  
  /*! 
    \fn void Final_State_Shower::InitTwojetTree(Tree *,double);
    \brief Initialise two toy-jets - mainly used for self-tests of
    the Final_State_Shower etc..
  */

  //-----------------------------------------------------------------------
  //------------------------ Evolution of the Shower ----------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn bool Final_State_Shower::EvolveJet(Tree *,Knot *);
    \brief Evolves the parton system radiated off an given knot

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



  /*!
    \fn bool Final_State_Shower::FillBranch(Tree *,Knot *,int);
    \brief  Try to branch a mother knot with given t and z.

    Try to branch a mother knot with given t and z :
    One of the daughter gets diced down more and the daughters system is
    checked for a kinematic fit. Then 1 is returned when a suitable branch 
    has been found, i.e. kinematically allowed daughters. Note, 
    that depending on the integer flag, either only the mothers z (i.e. the energy 
    fractions of the daughters) or the daughters momenta (if the daughters are from a 
    ME) get shuffled and eventually have to be restored.
  */

  //-----------------------------------------------------------------------
  //------------------- Service for the Branchings ------------------------
  //----------------------------------------------------------------------- 
  /*!
    \fn     bool Final_State_Shower::ChooseDaughter(Knot *);
    \brief New daughter to dice {0 or 1} for {d1 or d2}
      is chosen according to who's farer off-shell.
      Constraints on daughters apply if one of them
      is already inactive, i.e. its t is set on t_out.
  */


  /*!
    \fn void Final_State_Shower::InitDaughters(Tree *,Knot *,APHYTOOLS::Flavour *,bool);
    \brief
      Initializes new daughters. If knots do not exist so far,
      they are initialized including there mother-daughter relation
      and the flavour they carry. If they are already initialised,
      their status is reset to "active", and their t and E^2 are reset
  */
  

  /*!
    \fn   void Final_State_Shower::ResetDaughters(Knot *);
    \brief Calls to reset the left and the right knot of the mother
  */


  /*!
    \fn void Final_State_Shower::Reset(Knot *);
    \brief Deletes recursively the daughters and next generations of 
    each knot, partons are set inactive (leading to empty entries
    in the partonlist), they status of the knots is set to inactive,
    their t to zero and their daughters get deleted.
  */
  
  /*!
    \fn bool Final_State_Shower::ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2);
    \brief make sure number of jet from hard event is not reduced!!!
  */

  /*!
    \fn bool Final_State_Shower::ExtraJetCheck();
    \brief make sure number of jet from hard event is not reduced using ALL ini_partons
  */

}
#endif



