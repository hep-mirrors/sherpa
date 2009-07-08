#ifndef APACIC_Showers_Spacelike_Sudakov_H
#define APACIC_Showers_Spacelike_Sudakov_H

#include "AddOns/Apacic++/Showers/Backward_Splitting_Group.H"
#include "ATOOLS/Math/Random.H"

namespace ATOOLS { class Data_Reader; }

namespace PDF {
  class PDF_Base;
  class Remnant_Base;
}

namespace APACIC {
  class Sudakov_Tools;
  class Spacelike_Kinematics;
  class Knot;

  // dicing z and t for spacelike (initial state) shower evolution
  class Spacelike_Sudakov : public Backward_Splitting_Group {
    // global parameters
    ATOOLS::Mass_Selector *p_rms;
    Sudakov_Tools        * p_tools;
    Spacelike_Kinematics * p_kin;
    PDF::PDF_Base        * p_pdfa;
    PDF::Remnant_Base    * p_remnant;
    double m_pt2min, m_s_hadron, m_t0, m_emin;

    int    m_cpl_scheme,m_pdf_scheme, m_ordering_scheme;
    static size_t s_nem, s_maxem;

    double m_pdf_fac, m_pdf_scale_fac, m_cpl_scale_fac;
    // mother
    double m_E2, m_x, m_z, m_t, m_ta, m_tc, m_phi, m_pt2, m_rbmax;
    ATOOLS::Flavour m_inflav;

    std::vector<int> m_vetos;

    void Add(Splitting_Function * spl) ;
    double CrudeInt(double _zmin, double _zmax);    
    inline void SelectOne() { p_selected->SelectOne(); }

    void ProduceT();  
    bool Veto(Knot *); 
    bool CplVeto();
    bool OrderingVeto(Knot *);
    bool MassVeto();
    bool RemnantVeto(Knot *);
    void UniformPhi() { m_phi =  2.*M_PI*ATOOLS::ran.Get(); };
  public:
    Spacelike_Sudakov(PDF::PDF_Base *,Sudakov_Tools *, Spacelike_Kinematics *,
		      double,ATOOLS::Data_Reader *,const size_t); 
    ~Spacelike_Sudakov();
    bool Initialize();
    bool Dice(Knot * mother,double sprime,double m2p,int & extra_pdf);
    void AcceptBranch(const Knot *const mo);
    bool PTVeto(Knot *);

    static double Q02(ATOOLS::Flavour const & flav) { 
      return ATOOLS::sqr(flav.Mass()); 
    } // !!!
    static double MinE2(ATOOLS::Flavour const & flav) { 
      return ATOOLS::sqr(flav.Mass()); 
    } //!!!mass

    inline void SetRemnant(PDF::Remnant_Base *const remnant) { p_remnant=remnant; }

    inline void SetMaxScale(const double &scale) { m_s_hadron=scale; }
    inline void SetMinEnergy(const double &emin) { m_emin=emin;      }

    inline void SetPDFScaleFactor(const double &fac)      { m_pdf_scale_fac=fac; }
    inline void SetCouplingScaleFactor(const double &fac) { m_cpl_scale_fac=fac; }

    inline void SetCouplingScheme(const int &scheme) { m_cpl_scheme=scheme;      }
    inline void SetPDFScheme(const int &scheme)      { m_pdf_scheme=scheme;      }
    inline void SetOrderingScheme(const int &scheme) { m_ordering_scheme=scheme; }

    inline double PT2() const { return m_pt2; }

    inline void SetRBMax(const double &rb) { m_rbmax=rb; }

    inline double RBMax() const { return m_rbmax; }

    inline static void SetMaxEmissions(const size_t &n) { s_maxem=n; }
    inline static void SetEmissions(const size_t &n)    { s_nem=n;   }
    inline static void AddEmission()                    { ++s_nem;   }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_rms=ms; }

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

  /*!
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
    \fn static double Spacelike_Sudakov::Q02(ATOOLS::Flavour const & flav);
    \brief returns flavour dependent Q0
  */

  /*!
    \fn static double Spacelike_Sudakov::MinE2(ATOOLS::Flavour const & flav)
    \brief returns flavour dependent Emin
  */



} // namespace APACIC
#endif // Spacelike_Sudakov_H
