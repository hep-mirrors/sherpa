#ifndef APACIC_Showers_Initial_State_Shower_H
#define APACIC_Showers_Initial_State_Shower_H

#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "AddOns/Apacic++/Showers/Spacelike_Kinematics.H"
#include "AddOns/Apacic++/Showers/Spacelike_Sudakov.H"
#include "AddOns/Apacic++/Showers/Final_State_Shower.H"
#include "AddOns/Apacic++/Showers/Timelike_Kinematics.H"
#include "AddOns/Apacic++/Showers/Sudakov_Tools.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Apacic++/Main/Tree.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace APACIC {

  class Initial_State_Shower {
  private:

    Final_State_Shower    *p_fin;
    Sudakov_Tools         *p_tools;
    Spacelike_Kinematics  *p_kin;
    Spacelike_Sudakov    **p_suds;

    Tree **p_istrees, *p_fstree;
    Knot  *p_k1, *p_k2;

    double m_t0, m_sprime, m_pt2_1, m_pt2_2, m_th_1, m_th_2;
    int    m_allowed, m_extra_pdf[2], m_to_be_diced[2], m_tree1;

    ATOOLS::Vec4D    m_lab;
    ATOOLS::Poincare m_cmsrot, m_labboost;

    ATOOLS::Blob *m_bl_meps_is, *m_bl_meps_fs;

    void SetDirection(Knot *const k);
    void SetVetoScales(Knot *const d);

    void ChooseMother(int &ntree1,Knot *&k1,Knot *&k2,bool init=false);

    bool InitializeSystem(Tree **const trees,int tree1,Knot *k1,Knot *k2);
    void InitTwoTrees(Tree ** trees,double E2);
    int  EvolveSystem(Tree **const trees,int tree1,Knot *k3,Knot *k5);
    int  FillBranch(Tree **const trees,const int tree1,
		    Knot *const active,Knot *const partner);

    void FillMotherAndSister(Tree * tree,Knot * k,ATOOLS::Flavour * k_flavs);
    void SetStartingConditions(double k1,double k2,int flag);

  public:

    // constructor
    Initial_State_Shower(PDF::ISR_Handler *const isr, 
			 Final_State_Shower *const fin,
			 MODEL::Model_Base *const model,
			 ATOOLS::Data_Reader *const read);

    // destructor
    ~Initial_State_Shower();

    // member functions
    void SingleExtract(ATOOLS::Blob *jet,Knot *const knot);

    int  DoKinematics();

    int  PerformShower(Tree **const initrees,Tree *const fintree);

    ATOOLS::Vec4D GetMomentum(Knot *const mo,int &number);

    void OutputTree(Tree * tree);

    void SetMS(ATOOLS::Mass_Selector *const ms);

    static void SetColours(Knot * k);

    // inline functions
    inline double CutOff() { return m_t0; }

    inline ATOOLS::Vec4D GetLab() { return m_lab; }

    inline void InitShowerPT(const double &pt2max)
    {
      m_pt2_1=m_pt2_2=pt2max;
      m_th_1=m_th_2=M_PI;
    }

  };

  /*!
    \file 
    \brief contains the class APACIC::Initial_State_Shower
  */

  /*!
    \class Initial_State_Shower 
    \brief is the main class for simulation of initial state QED and/or QCD radiation
  */

  /*! 
    \fn    Initial_State_Shower::Initial_State_Shower(PDF::ISR_Handler *,Final_State_Shower *)
    \brief Standard constructor
  */

  /*!
    \fn    Initial_State_Shower::~Initial_State_Shower()
    \brief destructor, deletes sudakov, kinematics and sudakov helpers.
  */  

  /*!
    \fn    bool Initial_State_Shower::PerformShower(Tree **, int)
    \brief Performs the Initial State Shower.
    
    This routine is called in order to start the shower evolution, employing
    InitializeSystem().
    After the evoluton some consistency checks are done.

    \warning
    Make sure, you initialized the two trees and the starting condition beforehand
    with something like InitShowerPT !
  */

  /*!
    \fn    void  Initial_State_Shower::InitShowerPT(double pt2max)
    \brief To start the shower from outside. 
  
    Hand over the maximal pt of the first branches in the backward evolution - usually
    given by    pt2_max = (2 * s * t *u) / (s^2 + t^2 + u^2) 
    in 2->2 processes.
  */

  /*!
    \fn bool Initial_State_Shower::InitializeSystem(Tree **,Knot *,Knot *)
    \brief determines starting points for the shower

    Tries to initialize the incoming system recursively along the following steps: 
    First, if any of the two daughters has a mother, i.e. is an internal line of
    a ME, InitializeSystem is called with this mother as an argument.
    If two motherless daughters are found, their spacelike virtualities are selected
    via FillBranch(). Then, this two particle  system is brought to its own cm system 
    and oriented along the z-axis. The evolution is performed by EvolveSystem().
    
    It might happen that at some stage some timelike evolution does not work out,
    then InitializeSystem returns 0 and the system has to be re-initialitzed.
  */


  /*!
    \fn bool Initial_State_Shower::EvolveSystem(Tree **,Knot *,Knot *)
    \brief    Evolves the system recursively.

    In each step a system consisting of two knots of different incoming trees
    is evolved by going a step back. The parton with larger virtuality is supposed 
    to be closer to the hard interaction and is selected for the step back.
    Its status is set on passive (0, no more radiation) and energies are
    selected (the sisters energy may change, though, if it goes off-shell, i.e.
    initiates a timelike shower.). Then the branch with the previous knot and
    the remaining partner is filled via FillBranch and the final 
    state evolution of the sister is performed calling 
    Final_State_Shower::FirstTimelikeFromSpacelike(), respectively. 
    If the kinematics called via Spacelike_Kinematics::DoKinematics() does not 
    work out EvolveSystem returns a 0 and the evolution terminates to be restarted anew.
    Otherwise, i.e. in case the kinematics works out, sprime is rescaled and
    the recursion continues.      
  */
  
  /*!
    \fn bool Initial_State_Shower::FillBranch(Tree **,Knot *,Knot *,int)
    \brief determines a new virtuality for a given knot

    Here, a new spacelike \e t for a specific knot is selected via
    Spacelike_Sudakov::Dice(). If this \e t is larger than the infrared cut-off 
    a mother and sister are initialized by FillMotherAndSister() with
    Flavours obtained from the Spacelike_Sudakov. The starting scale for the timelike
    shower is calculated with help of method provided by Spacelike_Kinematics 
    and the sister's t is set accordingly. If the selected is smaller than the 
    infrared cut-off, the active knot is set passive (stat = 0) and its t is
    set to zero.

    In case, the maximal allowed m^2 for the timelike branch is negative, a 0 is returned.
  */

  /*!
    \fn void Initial_State_Shower::FillMotherAndSister(Tree *,Knot *,ATOOLS::Flavour *)
    \brief extends a tree by the given information

    For a Knot beloning to a given tree a mother and a sister are initalized
    with the flavours as given. The basic kinematics and the mutual relationships 
    are established and colurs are set.
  */

  /*!
    \fn void Initial_State_Shower::SetColours(Knot *)
    \brief Selects colours for Initial_State_Shower::FillMotherAndSister.
  */

  /*!
    \fn void Initial_State_Shower::SetStartingConditions(double,double,int) ;
    \brief Set starting conditions for the two incoming showers.

    The int flag selects wether the starting condition is formulated in terms
    of transverse momenta (flag = 0) or in terms of angles for the first 
    branching (flag = 1).
  */


  /*!
    \fn void Initial_State_Shower::ExtractPartons(Knot *, ATOOLS::Blob *, 
                         ATOOLS::Blob_List *, ATOOLS::Particle_List *, int)
    \brief fills blobs and parton list according to the result of the shower evolution
    
    After the shower is complete the outgoing knots of the tree are put into
    the parton list. Moreover, blobs are filled.
  */

  /*!
    \fn bool Initial_State_Shower::TestShower(Tree **) ;
    \brief  Test the Initial State Shower with two toy-jets.

    This will initialize a test event. It is used to do checks of the initial
    state shower. Process will be something like  
    \f$ q \bar q \to  \gamma^* \f$  at a fixed energy scale. 
  */

  /*!
    \fn void Initial_State_Shower::InitTwoTrees(Tree **,double);
    \brief Two incoming trees are initialized to be used to test the shower routines.
  */


  /*! 
    \fn Knot * Initial_State_Shower::GetInitiator(int no);
    \brief Goes to the most previous knot of tree number \e no, i.e. the knot
    assumed to be the initiator of the full inital state shower.
  */

  /*! 
    \fn Knot * Initial_State_Shower::GetInitiator(Tree *);
    \brief Goes to the most previous knot of the tree, i.e. the knot
    assumed to be the initiator of the full initial state shower.
  */

  /*!
    \fn ATOOLS::Vec4D Initial_State_Shower::GetLab();
    \brief retruns a vetor to the LAB system
  */

  /*!
    \fn ATOOLS::Vec4D Initial_State_Shower::GetMomentum(Knot * mo,int & number);
    \brief adds all daughter momenta. Used for check four momentum conservation.
  */

  /*!
    \fn void Initial_State_Shower::OutputTree(Tree *);
    \brief prints out tree structure 
  */

} // namespace APACIC

#endif






