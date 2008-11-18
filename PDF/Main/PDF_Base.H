#ifndef PDF_Base_H
#define PDF_Base_H

#include "CXXFLAGS.H"
#include "Flavour.H"

namespace ATOOLS {

  class Integration_Info;

}

namespace PDF {

  class PDF_Base {

  private:

    static struct Box {
      std::vector<PDF_Base*> v_pdfp;
      Box() : v_pdfp() {}
      ~Box();
      std::size_t TrueEntryNumber() const;
    } s_box;

  protected:

    std::vector<PDF_Base*> m_copies;

    std::vector<ATOOLS::Flavour> m_partons;    
    ATOOLS::Flavour              m_bunch;

    std::string m_type;    

    double m_exponent, m_rescale;
    double m_xmin, m_xmax, m_q2min, m_q2max;

    std::vector<std::pair<ATOOLS::Flavour,double> > m_extracted;

    double m_fac_scale_factor;

    void SingleExtract(const ATOOLS::Flavour flavour,const double x);
    void SingleReset();

  public:

    // constructor
    PDF_Base();

    // destructor
    virtual ~PDF_Base();

    // member functions
    virtual void AssignKeys(ATOOLS::Integration_Info *const info) = 0;

    virtual bool   Collinear(const double kp2) const;
    virtual double Cut(const std::string &type);
    
    virtual void   Calculate(double x,double Q2) = 0;
    virtual double GetXPDF(const ATOOLS::Flavour flavour) = 0;

    virtual void      Output() = 0;
    virtual PDF_Base *GetCopy() = 0;

    virtual PDF_Base *GetBasicPDF();

    // inline functions
    inline double XMin() const  { return m_xmin;  }
    inline double XMax() const  { return m_xmax;  }
    inline double Q2Min() const { return m_q2min; }
    inline double Q2Max() const { return m_q2max; }

    inline const std::vector<ATOOLS::Flavour> &Partons() const { return m_partons; }
    inline const ATOOLS::Flavour               Bunch() const   { return m_bunch;   }

    inline const std::string &Type() const { return m_type;     }  
    inline double         Exponent() const { return m_exponent; }

    inline PDF_Base *operator->() { return GetBasicPDF(); }

    inline void Extract(const ATOOLS::Flavour flavour,const double x)
    { GetBasicPDF()->SingleExtract(flavour,x); }
    inline void Reset()
    { GetBasicPDF()->SingleReset(); }

  };
  
  /*!
    \class PDF_Base
    \brief is the base class for parton distributions and structure functions.

    This class is the base class for both parton distributions and structure functions
    (e.g. of leptons). It hosts some information about both the incoming bunch particle 
    and the partons that constitute it. Furthermore, the corresponding weights for 
    given energy fractions and scales are calculated and made available through this class.
  */
  /*!
    \var std::vector<ATOOLS::Flavour> PDF_Base::m_partons    
    The vector listing the constituents. for protons this list has usually eleven members,
    one gluon and the five lightest quarks and anti-quarks. For leptons, due to the very nature
    of the structure function, only the lepton in question is included as legitimate parton.
  */
  /*!
    \var ATOOLS::Flavour PDF_Base::m_bunch
    The incoming particle. 
  */
  /*!
    \var std::string PDF_Base::m_type
    The type of the PDF_Base. It is given by PDF_flavourname.
  */
  /*!
    \fn virtual PDF_Base * PDF_Base::GetCopy()
    A method to initialize another PDF as exact copy of the current one.
    This is needed for the initial state shower of APACIC.
  */
  /*!
    \fn void PDF_Base::Calculate(const double, const double)
    This calculates all parton densities inside a hadron or the structure function
    depending on the energy fraction \f$x\f$ and the scale squared \f$Q^2\f$. It should
    be noted that the result will be given - as usual - as \f$x\f$ times the function. 
  */
  /*!
    \fn double PDF_Base::GetXPDF(const ATOOLS::Flavour &)
    This returns x times one specific pdf out of the set of all calculated ones - it checks 
    whether the flavour passed over as argument is included in the list of legitimate constituents, 
    PDF_Base::m_partons and selects the appropriate one.
  */
  /*!
    \fn void PDF_Base::Output()
    Some default output for formatted output in run summary files. We'll have to work on
    such an option.
    \todo Work on output, maybe printing the PDFs ?
  */
} 

#endif
