#ifndef Cluster_Partons_Base_H
#define Cluster_Partons_Base_H

#include "Combine_Table_Base.H"
#include "NLL_Sudakov.H"
#include "Running_AlphaS.H"
#include "Jet_Finder.H"
#include "Matrix_Element_Handler.H"
#include "XS_Selector.H"
#include "Amegic.H"
#include "Exception.H"
#include <map>


namespace SHERPA {
  typedef std::map<ATOOLS::Flavour,ATOOLS::Flavour> Flavour_Map;

  class Cluster_Partons_Base: public ATOOLS::Terminator_Object {
  protected:
    Matrix_Element_Handler * p_me;
    MODEL::Running_AlphaS  * p_runas;
    ATOOLS::Jet_Finder     * p_jf;
    NLL_Sudakov            * p_fssud, * p_issud;
    Combine_Table_Base     * p_ct, * p_combi;
    ATOOLS::Blob           * p_blob;

    std::vector<double> m_qmin;

    int     m_njet, m_hardjets, m_maxqcdjets, m_minqcdjets;
    int     m_maxjetnumber, m_cut, m_max_decayjets;
    double  m_weight, m_asweight, m_scale, m_ycut, m_qmin_i, m_qmin_f, m_q2_f[2];
    double  m_q2_fss, m_q2_iss, m_q2_isjet, m_q2_fsjet, m_q2_amegic, m_jv_pt2;
    int     m_colors[4][2];

    int     m_isrmode, m_isrshoweron, m_fsrshoweron;
    int     m_sud_mode, m_bp_mode, m_as_mode, m_muf_mode;
    double  m_kfac, m_is_as_factor, m_fs_as_factor, m_me_as_factor;

    int     m_hard_nqcd, m_hard_nqed, m_nstrong;
    double  m_counts, m_fails, m_q2_hard;

    EXTRAXS::XS_Base * p_xs;

    // some statistics:
    std::vector<long>    m_events;
    std::vector<double>  m_meweight_sum, m_weight_sum, m_asweight_sum, m_sweight_sum;
    std::vector<double>  m_weight_sum_sqr, m_asweight_sum_sqr, m_sweight_sum_sqr;

    void SetQMin(Combine_Table_Base *const ct,
		 Combine_Table_Base *const ref,size_t &decayjets,
		 size_t &hardjets,double &hardqmin);

    Leg   **CreateLegs(int &,const int,const bool=false);
    bool    FillLegs(Leg *, AMEGIC::Point *, int &, int);
    void    CreateTables(Leg **,const int,const double,const double);

    int     Set2Colours(const int,const int,ATOOLS::Vec4D *,ATOOLS::Flavour *);
    int     Set3Colours(const int,const int,ATOOLS::Vec4D *,ATOOLS::Flavour *);
    int     Set4Colours(const int,const int,ATOOLS::Vec4D *,ATOOLS::Flavour *);

    void    PrepareTerminate();
    void    WriteOutWeights();
    virtual Combine_Table_Base *CreateTable(ATOOLS::Jet_Finder *,ATOOLS::Vec4D *,
					    Combine_Table_Base *,const int,const int) = 0;
  public :
    Cluster_Partons_Base(Matrix_Element_Handler *,ATOOLS::Jet_Finder * jf,const int,const int,const int,const int) ;
    virtual ~Cluster_Partons_Base();
    
    bool ClusterConfiguration(ATOOLS::Blob *,double=1.,double=1.);
    int  SetColours(ATOOLS::Vec4D *,ATOOLS::Flavour *);
    int  SetDecayColours(ATOOLS::Vec4D *,ATOOLS::Flavour *,int,int);
    void FixJetvetoPt2(double &);
    void JetvetoPt2(double &, double &,double &,double &);
    
    virtual EXTRAXS::XS_Base * GetXS(EXTRAXS::XS_Group *,ATOOLS::Flavour *) = 0;
    virtual int                SetColours(EXTRAXS::XS_Base *,ATOOLS::Vec4D *,ATOOLS::Flavour *) = 0;
    virtual void               CalculateWeight(const double &meweight) = 0;
    void SetQMin(Combine_Table_Base *const ref=NULL);

    ATOOLS::Flavour            Flav(int);
    ATOOLS::Vec4D              Momentum(int); 
    int                        Colour(int,int);
    Combine_Table_Base *       GetCombineTable();
    Flavour_Map *              GetFlavourMap();
    // inline functions
    inline double Weight() const      { return m_weight;    }
    inline double ISJetScale() const  { return m_q2_isjet;  }
    inline double FSJetScale() const  { return m_q2_fsjet;  }
    inline double YCut() const        { return m_ycut;      }
    inline int    OrderStrong() const { return m_hard_nqcd; }
    inline int    OrderEWeak()  const { return m_hard_nqed; }
    inline int    InSwaped() const    { return p_me->InSwaped(); }

    inline double FSShowerScale() const { return m_q2_fss;   }
    inline double ISShowerScale() const { return m_q2_iss;   }

    inline double FactorisationScale(const size_t &i) const { return m_q2_f[i]; }
    
    inline void SetMaxQCDJets(const int jets) { m_maxqcdjets=jets; }
    inline void SetMinQCDJets(const int jets) { m_minqcdjets=jets; }
  };

  /*!
    \file Cluster_Partons.H
    \brief Contains the class SHERPA::Cluster_Partons
  */

  /*!
    \class Cluster_Partons
    \brief main steering class for merging PS & ME
    
    This Class controls the merging of Parton Shower (PS) and
    Matrix Elements (ME). The main steps are the following:
      - ClusterConfiguration() : determines for a given 
        \f$2 \to n\f$ ME a parton shower history and a hard
        \f$2 \to 2\f$ Process
      - CalculateWeight() : evaluates the NLL Sudakov weight
        used for reweighting the ME kinematics
      - FillTrees() : determines the starting condition of the 
        PS, filling all information in the corresponding binary tree 
        structure
  */

  /*!
    \fn Cluster_Partons::Cluster_Partons(Matrix_Element_Handler * me, ATOOLS::Jet_Finder * jf, 
		    int maxjetnumber, int isrmode, int isrshoweron, int fsrshoweron)
    \brief Constructor, performs the basic initialization.

    The information needed for the initialization are:
     - the used Matrix_Element_Handler, to provide information like the included 
       Feynman diagrams
     - the used jet measure (class Jet_Finder)
     - the maximal possble number of outgoing particles, and
     - information whether the parton showers are available
  */

  /*!
    \fn  bool Cluster_Partons::ClusterConfiguration(ATOOLS::Blob * blob, double x1 ,double x2)
    \brief Applys a jet cluster algorithm one a given process until
           \f$2 \to 2\f$ Process is obtained

    The aim of this routine is to determine the best clustering in order to initialise the PS
    appropriately. A jet measure is employed to determine the pair of particles to be 
    clustered first. In order to disable unphysical parton histories, only those
    clusterings are possible, which have a corresponding Feynman diagram (in the 
    matrix element). 
    The jet measure to be used is set up during the initialization
    process. The default jet
    measure is \f$k_\perp\f$--cluster algorithm (cf. class ATOOLS::JetFinder).
    
    The clustering is performed with the help of the class Combine_Table.
  */

  /*!
    \fn  void Cluster_Partons::CalculateWeight(double hardscale, double asscale, double jetscale,double qm_i,double qm_f)
    \brief  Evaluates the NLL Sudakov factors

    Evaluates the sudakov weight for matching ME to PS. The Matrix Elements
    are reweighted with this weight.

    The sequence of clustering and the corresponding scales \f$q_i\f$ have been 
    evaluated ClusterConfiguration() with the help of the class Combine_Table. Here 
    this information is used to construct the sudakov weight, which 
    can be used to accept or reject events appropriately.

    The starting scale for the sudakov weight \e hardscale and the scale for 
    \f$\alpha_S\f$ factors \e asscale, was (in most cases) determined with call
    to SetColours(). The resolution scale \e jetscale is fixed by the used 
    jet measure used.

    The general form of the weight for an internal line is given by 
    \f[
       w = \frac{\Delta(q_i)}{\Delta(q_j)} \frac{\alpha_S(q_j)}{\alpha_S(q_0)} \,,
    \f]
    where the Sudakov form factors \f$\Delta(q_i)\f$ are calculated with the help of 
    the class NLL_Sudakov.

    Note that matrix element events with a maximal number of outgoing partices 
    have a reduced weight.
  */

  /*!
    \fn  void Cluster_Partons::FillTrees(APACIC::Tree ** ini_trees,APACIC::Tree * fin_tree)
    \brief creates tree structures to initialize the initial and final
           state parton showers

    Translate cluster sequence (as stored in Combine_Table's) into
    initial or final state shower seeds.

    Employs EstablishRelations(), DetermineColourAngles(), Point2Knot()
  */

  /*!
    \fn  int  Cluster_Partons::SetColours(ATOOLS::Vec4D * p, ATOOLS::Flavour * fl)
    \brief Determines a Colour Structure for a given hard \f$2 \to 2\f$ process
    and fixes a hard scale. 

    Determines a Colour Structure for a given hard \f$2 \to 2\f$ process
    and fixes a hard scale \e m_scale, and possibly a scale used in \f$\alpha_S\f$ factor
    \e m_asscale. In many cases \e m_scale and \e m_asscale are identical, only
    when three out of four particles are strong a difference might arise.

    This information can be used as initial condition for 
    CalculateWeight(). 
  
   \note This routine is only called, when the hard \f$2 \to 2\f$ process
    is not in available in the fast cross section library (cf. module EXTRAXS).
  */

  /*!
    \fn void Cluster_Partons::DetermineColourAngles(const std::vector<APACIC::Knot *> & knots)
    \brief determines the maximum angle between color connected partons of
    a hard \f$2 \to 2\f$ process.
    
    Determines angles to be used in angular ordering. 
    The colour angle is determined for initial state particles in the lab frame
    via boost along the z axis, whereas starting angles for finial state particles
    are determined in the current cms frame.
    The starting angles are stored in "thcrit" of each knot.
    Finially all momenta are restored.

    Utilizes ColourAngle().
  */

  /*!
    \fn void Cluster_Partons::EstablishRelations(APACIC::Knot * mo,APACIC::Knot * d1,APACIC::Knot * d2,int mode)

    Three cases (\e mode):
      -# two incoming partons from the hard \f$2 \to 2\f$ process
      -# two outgoing partons from a common mother 
      -# one incoming parton its mother and its sister
     .
     
     For the initialization of the incoming particles of 
     the hard \f$2 \to 2\f$ process the energy fractions \f$x_1\f$ and \f$x_2\f$
     are filled from information in the Combine_Table
     The starting scale of the shower, is given by the 
     virtual mass of the mother (since we have a shower evolution in virtualities).

     The two final state partices are initialized using
      APACIC::Final_State_Shower::EstablishRelations(mo,d1,d2) meaning:
     the more energetic parton is initialized with the angle and virtuality
     of the mother 
     the less energetic parton is initialized with angle and virtuality of
     the current branch.

     The one incoming parton its mother and its sister are initialized using
      APACIC::Initial_State_Shower::SetColours().
     Note, angle conditions will be fixed during the shower evolution.
  */

  /*!
    \fn double  Cluster_Partons::Weight()
    \brief returns the weight calculated in CalculateWeight()
  */
  /*!
    \fn double  Cluster_Partons::Scale()
    \brief returns the hardest scale (of the core process) determined in SetColours()
  */
  /*!
    \fn double  Cluster_Partons::AsScale()
    \brief return the scale associated with strong coupling (in the core process) as determined in SetColours()
  */
  /*!
    \fn ATOOLS::Flavour Cluster_Partons::Flav(int)
    \brief provides a simple access routine to the flavours of the core \f$2 \to 2\f$ process
  */
  /*!
    \fn ATOOLS::Vec4D   Cluster_Partons::Momentum(int)
    \brief provides a simple access routine to the momenta of the core \f$2 \to 2\f$ process
  */
  /*!
    \fn APACIC::Knot  * Cluster_Partons::Point2Knot(APACIC::Tree *,const Leg &,const ATOOLS::Vec4D &,char)
    \brief helper to translate information from Combine_Table to APACIC::Tree, used in FillTrees()
  */
  /*!
    \fn void    Cluster_Partons::CreateFlavourMap()
    \brief determines information to translate Feynman diagrams into shower histories

    When simulating many processes at once it can happen, that the same matrix element
    can be used to describe more than one process. In most cases it is due to a 
    simple change of the flavour of a massless fermion. In order to save memory and 
    computing time those processes are linked to "partner" processes.

    In order to translate the Feynman diagram of a "partner" processes, to feynman diagrams
    of the process in question a translation table (a map) is created.
  */

}
#endif
