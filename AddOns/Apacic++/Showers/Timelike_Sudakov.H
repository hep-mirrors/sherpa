#ifndef APACIC_Showers_Timelike_Sudakov_H
#define APACIC_Showers_Timelike_Sudakov_H

#include "AddOns/Apacic++/Showers/Splitting_Group.H"
#include "MODEL/Main/Model_Base.H"

namespace APACIC {

  class Knot;
  class Sudakov_Tools;
  class Timelike_Kinematics;

  class Timelike_Sudakov : public Splitting_Group {
  private:

    Sudakov_Tools       *p_tools;
    Timelike_Kinematics *p_kin;

    ATOOLS::Mass_Selector *p_rms;

    double m_pt2min, m_pt2max, m_t0; 
    double m_alphasmax, m_alphasfix; 
    double m_pt2min_qed, m_t0_qed, m_alpha_qed;

    int m_cpl_scheme, m_evolution_scheme, m_ordering_scheme;
    int m_mass_scheme, m_width_scheme, m_zrange_scheme;
    int m_mecorr_scheme,m_qed_mecorr_scheme,m_direct_photons,m_azimuthal_correlation;
    int m_last_veto, m_mode, m_shower, m_kt_scheme;
  
    double m_E2, m_t, m_tb, m_tc, m_z, m_phi, m_pt2;
    double m_maxpt, m_rscalefac, m_oldt, m_rbmax, m_fcs;

    size_t m_nem, m_maxem;

    const Knot *p_lemkn;

    ATOOLS::Simple_Polarisation_Info m_pol_b, m_pol_c;

    ATOOLS::Flavour m_inflav;

    double CrudeInt(double zmin,double zmax) ;

    void Add(Splitting_Function *spl);
    void SelectOne();

    void ProduceT(double ta, double m2);  

    bool Veto(Knot *const mo, double ta); 

    bool CplVeto(double t);
    bool OrderingVeto(Knot *const mo,double ta,double wa,double z);
    bool MEVeto(Knot *const mo,double t);
    bool MassVeto(double ta,double wa,double z);

    double UniformPhi() const;

    bool DiceT(Knot *const mother,Knot *const granny);

#ifdef CHECK_SPLITTINGS
    void CheckSplittings();
#endif

  public:

    // constructor
    Timelike_Sudakov(Timelike_Kinematics *const kin,
		     MODEL::Model_Base *const model);

    // destructor
    ~Timelike_Sudakov();

    // member functions
    void Init(const double fmed=0.0);
    bool Dice(Knot *const mother,Knot *const granny=NULL);
    void AcceptBranch(const Knot *const mo);
    void AcceptEmission(const Knot *const mo);

    bool MEWeight(Knot *const me);

    const ATOOLS::Simple_Polarisation_Info &GetPolB();
    const ATOOLS::Simple_Polarisation_Info &GetPolC();

    void SetPT2Min(const double &pt2);

    // inline functions
    inline void SetCouplingScheme(const int scheme)     
    { m_cpl_scheme=scheme; }
    inline void SetEvolutionScheme(const int scheme)     
    { m_evolution_scheme=scheme; }
    inline void SetOrderingScheme(const int scheme)     
    { m_ordering_scheme=scheme; }
    inline void SetMassScheme(const int scheme)         
    { m_mass_scheme=scheme; }
    inline void SetWidthScheme(const int scheme)        
    { m_width_scheme=scheme; }
    inline void SetMECorrectionScheme(const int scheme) 
    { m_mecorr_scheme=scheme; }
    inline void SetQEDMECorrectionScheme(const int scheme) 
    { m_qed_mecorr_scheme=scheme; }
    inline void SetQEDScheme(const int scheme)          
    { m_direct_photons=scheme; }
    inline void SetCorrelationScheme(const int scheme)  
    { m_azimuthal_correlation=scheme; }
    inline void SetKTScheme(const int scheme)  
    { m_kt_scheme=scheme; }

    inline void SetScaleFactor(const double &fac) { m_rscalefac=fac;   }
    inline void SetPT2MinQED(const double &pt2)   
    { m_t0_qed=4.0*(m_pt2min_qed=pt2); }

    inline int OrderingScheme() const { return m_ordering_scheme; }

    inline double ScaleFactor() const { return m_rscalefac; }
    inline double PT2Min() const      { return m_pt2min;    }

    inline void ResetN() { m_nem=0; }

    inline void SetNMax(const size_t &n) { m_maxem=n; }

    inline size_t GetN()    { return m_nem;   }
    inline size_t GetNMax() { return m_maxem; }

    inline const Knot *LastEmission() const { return p_lemkn; }

    inline void SetRBMax(const double &rb) { m_rbmax=rb; }
    inline void SetFCS(const double &fcs)  { m_fcs=fcs;  }

    inline double RBMax() const { return m_rbmax; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_rms=ms; }

  };// end of class Timelike_Sudakov

  /*!
    \file
    \brief contains the class APACIC::Timelike_Sudakov
  */

  /*!
    \class Timelike_Sudakov
    \brief Class for dicing z and t for timelike (final state) shower evolution
  */

  /*! 
    \fn Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics *);
    \brief Standard Constructor -

      loads of stuff to be removed ....
  */
  /*!
    \fn  Timelike_Sudakov::~Timelike_Sudakov();
    \brief Destructor
  */

  /*! 
    \fn bool Timelike_Sudakov::Dice(Knot * mother,Knot * granny=0);
    \brief Tries to generate t and z for a mother ( returns "1" if successful),
       needs virtuality of granny and energy and flavour of mother
  */
  
  /*! 
    \fn    void Timelike_Sudakov::Add(Splitting_Function * spl);
    \brief Initialises lists of splitting functions grouped by the flavour of the
           incoming leg A.
  */
  
  /*!
     \fn     double Timelike_Sudakov::CrudeInt(double _zmin, double _zmax) ;
     \brief Yields the crude integral of a group.
  */

     
  /*!
    \fn void Timelike_Sudakov::SelectOne();
    \brief Selects one specific mode for the splitting (to be called after CrudeInt)
  */

  /*! 
    \fn void Timelike_Sudakov::ProduceT();  
    \brief Starting with a t a new, a lower t will be produced and set.
    For quarks there's also the possibility to have q->q+photon. 
  */
  /*!
    \fn bool Timelike_Sudakov::Veto(Knot *);  
    \brief Organizes the vetoes to cure overshoot in approximate Sudakov formfactor
    used in ProduceT, returns a 0 or 1 if no veto or a veto is thrown.
    Calls the specific vetos after some trivial tests.
  */
  
  /*!
    \fn     bool  Timelike_Sudakov::CplVeto();
    \brief Corrects for overshooting in taking a crude alpha_S
  */

  /*!
    \fn     bool  Timelike_Sudakov::AngleVeto(Knot *);
    \brief Veto from imposed angular ordering on subsequent emissions
  */
  
  /*!
    \fn   bool Timelike_Sudakov::MEVeto(Knot *);
    \brief applies ME corrections to first correction

    Depending on ordering schemes,emissions get corrected on
    the ME kinematics of e+e- -> qqg
    if "first" branch perform ME - Correction for gluon radiation
    \verbatim
             (t',z') *      
               /     \      x_i = 2 E_i / sqrt(t')
            (t,z) *   \
             / \       \
            1   3       2
    \endverbatim
  */
  
  /*! 
    \fn bool Timelike_Sudakov::MassVeto();
    \brief checks for \e z being constient with physical range

    First veto from condition on cos(theta_bc) to be physical, condition:
    \f$ z_min < z < z_max\f$,   \f$z_min,max = 1/2 (1 -+ p/E)\f$
    Here, there's a constrained and an unconstrained choice.
    Second veto from modified splitting function for heavy quark a la
    Catani. 
  */
  
  /*!
    \fn bool Timelike_Sudakov::JetVeto(Knot *);  
    \brief considers only branching that do not constitute an independen jet
  */
    

  /*!
    \fn     void   Timelike_Sudakov::UniformPhi()
    \brief Uniform distribution of phi, angular correlations to be imposed.
  */


  /*!
    \fn void Timelike_Sudakov::CheckSplittings()
    \brief Perfroming consistency checks; for Internal use only
  */

}// end of namespace APACIC

#endif
