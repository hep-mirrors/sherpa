#ifndef Cluster_Partons_H
#define Cluster_Partons_H


#include "Combine_Table.H"
#include "NLL_Sudakov.H"

#include "Tree.H"
#include "Knot.H"

#include "Process_Base.H"
#include "XS_Base.H"
#include "Jet_Finder.H"

#include "Blob.H"

namespace SHERPA {
  class Matrix_Element_Handler;

  class Cluster_Partons {
    Combine_Table  * p_ct;
    Combine_Table  * p_combi;

    NLL_Sudakov    * p_sud;

    //    AMEGIC::Process_Base   * p_lastproc;
    Matrix_Element_Handler  * p_me;
    ATOOLS::Jet_Finder   * p_jf;
    ATOOLS::Blob         * p_blob;

    int     m_maxjetnumber;
    double  m_weight;

    double  m_scale;
    int     m_colors[4][2];

    int     m_sud_mode;
    int     m_isron, m_fsron;

    // some statistics:
    long    * p_events;
    double  * p_weight_sum;
    double  * p_weight_sum_sqr;


    bool            FillLegs(Leg *, AMEGIC::Point *, int &, int);
    APACIC::Knot  * Point2Knot(APACIC::Tree *,const Leg &,const ATOOLS::Vec4D &,char);
    ATOOLS::Vec4D Momentum(APACIC::Knot *, int &);

    bool    IsColourConnected(ATOOLS::Parton *,ATOOLS::Parton *);
    double  ColourAngle(const std::vector<APACIC::Knot *> &, const int);
    void    EstablishRelations(APACIC::Knot *,APACIC::Knot *,APACIC::Knot *,int);
    void    DetermineColourAngles(const std::vector<APACIC::Knot *> &);

  public :
    Cluster_Partons(Matrix_Element_Handler * me, ATOOLS::Jet_Finder * jf, 
		    int maxjetnumber, int isron, int fsron) ;
    ~Cluster_Partons();
    bool ClusterConfiguration(ATOOLS::Blob *,double,double);
    void CalculateWeight(double,double);
    void FillTrees(APACIC::Tree **,APACIC::Tree *,EXTRAXS::XS_Base *);

    int     SetColours(ATOOLS::Vec4D *, ATOOLS::Flavour *);
    double  Weight();
    double  Scale();

    ATOOLS::Flavour Flav(int);
    ATOOLS::Vec4D    Momentum(int); 
  };

  inline double  Cluster_Partons::Weight() { return m_weight; }
  inline double  Cluster_Partons::Scale()  { return m_scale; }


  /*!
    \class Cluster_Partons
    \brief main steering class for merging PS & ME
    
    This Class controls the merging of Parton Shower (PS) and
    Matrix Elements (ME). The main steps are the following:
      - ClusterConfiguration() : determines for a given 
        2 \to n ME a parton shower history and a hard
        2 \to 2 Process
      - CalculateWeight() : evaluates the a NLL Sudakov weight
        used for reweighting the ME kinematics
      - FillTrees() : determines the starting condition of the 
        PS, filling all information in the corresponding binary tree 
        structure
  */

  /*!
    \fn  bool Cluster_Partons::Cluster_Configuration(ATOOLS::Blob *)
    \brief Applys a jet cluster algorithm one a given process until
           2 \to 2 Process is obtained

    Determine best clustering, in order to initialise the PS there for
    a jet algorithm (given at initialization time is used). The actual
    work is done in Combine_Table.

  */

  /*!
    \fn  void Cluster_Partons::CalculateWeight(double hardscale ,double jetscale);
    \brief  Evaluates the NLL Sudakov factors

    Evaluates the sudakov weight for matching ME to PS. The Matrix Elements
    are reweighted with this weight.

  */

  /*!
    \fn  void Cluster_Partons::FillTrees(APACIC::Tree **,APACIC::Tree *,EXTRAXS::XS_Base *)
    \brief creates tree structures to initialize the parton initial and final
           state parton showers

    Translate cluster sequence (as stored in Combine_Table's) into
    initial or final state shower seeds.
    Employs EstablishRelations(), DetermineColourAngles()
  */

  /*!
    \fn  int  Cluster_Partons::SetColours(ATOOLS::Vec4D *, ATOOLS::Flavour *)
    \brief Determines a Colour Structure for a given hard 2 \to 2 process
    and fixes a hard scale. These can be used as initial condition for 
    CalculateWeight()
  */

  /*!
    \fn void DetermineColourAngles(const std::vector<APACIC::Knot *> &)
    
    First step boost all knots in thier cms frame, rotates incoming on z-axis, and 
    determines angles to be used in angular ordering utilizing ColourAngle():
    The colour angle is determined for initial state particles in the lab frame
    via boost along the z axis, whereas starting angles for finial state particles
    are determined in the current cms frame.
    The starting angles are stored in "thcrit" of each knot.
    Finially all momenta are restored.
  */

  /*!
    \fn void EstablishRelations(APACIC::Knot *,APACIC::Knot *,APACIC::Knot *,int)

    Three cases:
      -# two incoming partons from the hard 2 \to 2 process
      -# two outgoing partons from a common mother
      -# one incoming parton its mother and its sister
     .
     
     For the initialization of the incoming particles of 
     the hard 2 \to 2 process the energy fractions x_1 and x_2
     are filled from information in the Combine_Table
     The stating scale of the shower, is given by the 
     virtual mass of the mother (since we have a shower in virtualities).

     The two final state partices are initialized using
      APACIC::Final_State_Shower::EstablishRelations(mo,d1,d2) meaning:
     the more energetic parton is initialized with the angle and virtuality
     of the mother 
     the less energetic parton is initialized with angle and virtuality of
     the current branch

     The one incoming parton its mother and its sister are initialized using
      APACIC::Initial_State_Shower::SetColours
     Note, angle conditions will be fixed during the shower evolution!


       
  */


}
#endif
