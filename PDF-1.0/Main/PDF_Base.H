#ifndef PDF_Base_H
#define PDF_Base_H

#include "Flavour.H"

namespace PDF {
  class PDF_Base {
  protected :
    std::vector<APHYTOOLS::Flavour> m_partons;    
    APHYTOOLS::Flavour m_bunch;
    std::string        m_type;    
  public:
    virtual ~PDF_Base() {}
    virtual void Calculate(const double x, const double Q2) = 0;
    virtual double GetXPDF(const APHYTOOLS::Flavour &)      = 0;

    virtual void                            Output()        = 0;
    virtual PDF_Base                      * GetCopy()       = 0;
    virtual double                          GetXMin()       = 0;
    virtual double                          GetXMax()       = 0;
    virtual double                          GetQ2Min()      = 0;
    virtual double                          GetQ2Max()      = 0;
    virtual APHYTOOLS::Flavour              Bunch()         = 0;
    virtual std::vector<APHYTOOLS::Flavour> Partons()       = 0;
    virtual double                          Exponent()     { return 1.; }
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
    \var std::vector<APHYTOOLS::Flavour> PDF_Base::m_partons    
    The vector listing the constituents. for protons this list has usually eleven members,
    one gluon and the five lightest quarks and anti-quarks. For leptons, due to the very nature
    of the structure function, only the lepton in question is included as legitimate parton.
  */
  /*!
    \var APHYTOOLS::Flavour PDF_Base::m_bunch
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
    \fn double PDF_Base::GetXPDF(const APHYTOOLS::Flavour &)
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
