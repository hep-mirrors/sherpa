#ifndef Initial_State_Shower_H
#define Initial_State_Shower_H

#include "ISR_Handler.H"
#include "Spacelike_Kinematics.H"
#include "Spacelike_Sudakov.H"
#include "Final_State_Shower.H"
#include "Timelike_Kinematics.H"
#include "Sudakov_Tools.H"
#include "Parton_List.H"
#include "Tree.H"
#include "Run_Parameter.H"

namespace APACIC {
  class Initial_State_Shower {
  private:
    Spacelike_Kinematics *  kin;
    Spacelike_Sudakov    ** suds;
    Final_State_Shower   *  fin;
    Sudakov_Tools        *  tools;

    bool   jetveto;
    double sprime;
    double t0;

    double pt2_1,pt2_2,th_1,th_2;

    int    allowed;

    AMATOOLS::Vec4D lab;

  public:
    Initial_State_Shower(ISR::ISR_Handler *, Final_State_Shower *);// for testing purposes
    ~Initial_State_Shower();
    bool PerformShower(Tree **,bool);
    void InitShowerPT(double pt2max);
    void ExtractPartons(Knot *,APHYTOOLS::Blob *,
			APHYTOOLS::Blob_List *,APHYTOOLS::Parton_List *,int);
    bool TestShower(Tree **) ;

    AMATOOLS::Vec4D GetLab();
    AMATOOLS::Vec4D GetMomentum(Knot * mo,int & number);

    Knot * GetInitiator(int);
    Knot * GetInitiator(Tree *);

  private:
    bool InitializeSystem(Tree **,Knot *,Knot *);
    void InitTwoTrees(Tree **,double);
    bool EvolveSystem(Tree **,Knot *,Knot *);
    bool FillBranch(Tree **,Knot *,Knot *,int);

    void FillMotherAndSister(Tree *,Knot *,APHYTOOLS::Flavour *);
    void SetColours(Knot *);
    void SetStartingConditions(double,double,int) ;

    void OutputTree(Tree *);
  };

  //----------------------------------------------------------------------
  //    inline functions
  //----------------------------------------------------------------------


  inline void Initial_State_Shower::SetStartingConditions(double k1,double k2,int flag) 
  {
    switch (flag) {
    case  1 :  th_1  = k1; th_2  = k2; break;
    default :  pt2_1 = k1; pt2_2 = k2; break;
    };
  }

  inline AMATOOLS::Vec4D Initial_State_Shower::GetLab() { return lab; }

  inline Knot * Initial_State_Shower::GetInitiator(Tree * tree) 
  {
    return tree->GetInitiator();
  }

  inline AMATOOLS::Vec4D Initial_State_Shower::GetMomentum(Knot * mo,int & number) 
  {
    if (mo->left) return GetMomentum(mo->left,number) + GetMomentum(mo->right,number);
    number++;
    return mo->part->Momentum();
  }
} // namespace APACIC








/*!
  \class Initial_State_Shower 
*/

/*! 
  \fn    Initial_State_Shower(ISR::ISR_Handler *,Final_State_Shower *)
  \brief Standard constructor
*/
  
/*!
  \fn    bool PerformShower(Tree **);
  \brief Performs the Initial State Shower.

  Make sure, you initialized the two trees and the starting condition beforehand
  with something like InitShowerPT !
*/

/*!
  \fn    void  InitShowerPT(double pt2max);
  \brief To start the shower from outside. 
  
  Hand over pointers to the two trees to be filled
  and the maximal pt of the first branches in the backward evolution - usually
  given by    pt2_max = (2 * s * t *u) / (s^2 + t^2 + u^2) 
  in 2->2 processes.
*/

/*!
  \fn void ExtractPartons(Knot *,Blob *,Blob_List *,Parton_List *,int);
  After the shower is complete the outgoing knots of the tree are put into
  the parton list. Moreover, blobs are filled.
*/

/*!
  \fn bool TestShower(Tree **) ;
  Test the Initial State Shower with two toy-jets.
  Process will be something like   q bar q -> gamma*     at a fixed
  energy scale. 
*/

/*! 
  \fn Knot * GetInitiator(int);
  Goes to the most previous knot of tree number no, i.e. the knot
  assumed to be the initiator of the full inital state shower.
*/

/*! 
  \fn Knot * GetInitiator(Tree *);
  Goes to the most previous knot of the tree, i.e. the knot
  assumed to be the initiator of the full inital state shower.
*/

/*!
  \fn bool InitializeSystem(Tree **,Knot *,Knot *);
  Tries to initialize the incoming system recursively along the following steps: 
  First, if any of the two daughters has a mother, i.e. is an internal line of
  a ME, InitializeSystem is called with this mother as an argument.
  If two motherless daughters are found, their spacelike virtualities are selected
  via Initial_State_Shower::FillBranch. Then, this two particle 
  system is brought to its own cm system and oriented along the z-axis.
  The evolution is performed by Initial_State_Shower::EvolveSystem.
  It might happen that at some stage some timelike evolution does not work out,
  then InitializeSystem returns 0 and the system has to be re-initialitzed.
*/

/*!
  \fn void InitTwoTrees(Tree **,double);
  Two incoming trees are initialized to be used to test the shower routines.
*/

/*!
  \fn bool EvolveSystem(Tree **,Knot *,Knot *);
  Evolves the system recursively.
  In each step a system consisting of two knots of different incoming trees
  is evolved by going a step back. The parton with larger virtuality is supposed 
  to be closer to the hard interaction and is selected for the step back.
  Its status is set on passive (0, no more radiation) and energies are
  selected (the sisters energy may change, though, if it goes off-shell, i.e.
  initiates a timelike shower.). Then the branch with the previous knot and
  the remaining partner is filled via Initial_State_Shower::FillBranch and 
  the final state evolution of the sister is performed calling 
  Final_State_Shower::FirstTimelikeFromSpacelike, respectively. 
  If the kinematics called via Spacelike_Kinematics::DoKinematics does not 
  work out EvolveJet returns a 0 and the evolution terminates to be restarted anew.
  Otherwise, i.e. in case the kinematics works out, sprime is rescaled and
  the recursion continues.      
*/

/*!
  \fn bool FillBranch(Tree **,Knot *,Knot *,int);
  Here, a new spacelike t for a specific knot is selected via
  Spacelike_Sudakov::Dice. If this t is larger than the infrared cut-off a mother 
  and sister are initialized by Initial_State_Shower::FillMotherAndSister with
  Flavours obtained from the Sudakov. The starting scale for the timelike
  shower is calculated with help of method provided by Spacelike_Kinematics 
  and the sister's t is set accordingly. If the selected is smaller than the 
  infrared cut-off, the active knot is set passive (stat = 0) and its t is
  set to zero.
  In case, the maximal allowed m^2 for the timelike branch is negative, a 0 is returned.
*/

/*!
  \fn void FillMotherAndSister(Tree *,Knot *,APHYTOOLS::Flavour *);
  For a Knot beloning to a given tree (int) a mother and a sister are initalized
  with the flavours as given. The kinematics and the mutual relationships are established 
  and colurs are set.
*/

/*!
  \fn void SetColours(Knot *);
  Selects colours for Initial_State_Shower::FillMotherAndSister.
*/

/*!
  \fn void SetStartingConditions(double,double,int) ;
  Set starting conditions for the two incoming showers.
  The int flag selects wether the starting condition is formulated in terms
  of transverse momenta (flag = 0) or in terms of angles for the first 
  branching (flag = 1).
*/


#endif // Initial_State_Shower.H






