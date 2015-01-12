#ifndef PDF_MRST_PDF_MRST99_H
#define PDF_MRST_PDF_MRST99_H

#include <iostream>
#include "PDF/Main/PDF_Base.H"
#include "PDF/MRST/mrst99.h"

namespace PDF {

  class PDF_MRST99 : public PDF_Base {
  private:

    static c_mrst * p_proton;
    s_partoncontent m_content;
    std::string     m_path;
    int             m_anti;
    bool            m_overscaled;

  public:

    PDF_MRST99(const ATOOLS::Flavour,const int); 

    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour); 

  };

  /*!
    \class PDF_MRST99
    \brief An interface to the C++-version of the MRST99 sets.
    
    This is an interface of the SHERPA package to the 
    <A HREF="http://durpdg.dur.ac.uk/hepdata/mrstc++.html>C++-version of the MRST99<\A>
    sets of parton distribution functions by Jeppe Andersen. 
  */
  /*!
    \var static c_mrst * PDF_MRST99::proton
    A pointer to the C++ object of the MRST99 set.
  */
  /*!
    \var s_partoncontent PDF_MRST99::m_content
    Entries for the various PDFs.
  */
  /*!
    \var std::string PDF_MRST99::m_path
    The path for the grid to interpolate.
  */
  /*!
    \var int PDF_MRST99::m_set
    The set number of MRST99. In this case the sets denote:
    m_set   grid   comment               L(4)/MeV  a_s(M_Z)  grid#1
    - 1     COR01  central gluon, a_s    300       0.1175    0.00524 
    - 2     COR02  higher gluon          300       0.1175    0.00497 
    - 3     COR03  lower gluon           300       0.1175    0.00398  
    - 4     COR04  lower a_s             229       0.1125    0.00585  
    - 5     COR05  higher a_s            383       0.1225    0.00384  
    - 6     COR06  quarks up             303.3     0.1178    0.00497  
    - 7     COR07  quarks down           290.3     0.1171    0.00593  
    - 8     COR08  strange up            300       0.1175    0.00524  
    - 9     COR09  strange down          300       0.1175    0.00524  
    - 10    C0R10  charm up              300       0.1175    0.00525  
    - 11    COR11  charm down            300       0.1175    0.00524  
    - 12    COR12  larger d/u            300       0.1175    0.00515  
  */
  /*!
    \var int PDF_MRST99::m_anti
    To denote particles and antiparticles, m_anti = 1, -1, respectively.
  */
  /*!
    \fn PDF_MRST99::PDF_MRST99(const ATOOLS::Flavour,const int,const std::string); 
    In the constructor the mrst99 class is initialized and the list of partons is filled.
  */
  /*!
    \fn void PDF_MRST99::Calculate(const double, const double)
    This calculates all parton densities inside the proton depending on the energy fraction 
    \f$x\f$ and the scale squared \f$Q^2\f$. It should be noted that the result will be given - 
    as usual - as \f$x\f$ times the function. The calculation consists of calling the mrst99
    method of the c_mrst class and by updating the parton densities of m_content accordingly. 
  */
  /*!
    \fn double PDF_Base::GetXPDF(const ATOOLS::Flavour &)
    This returns x times one specific pdf out of the set of all calculated ones - it checks 
    whether the flavour passed over as argument is included in the list of legitimate constituents, 
    PDF_Base::m_partons and selects the appropriate one. In this case, the corresponding parton
    density might consist of sea and valence contributions or of sea contributions alone.
  */
  /*!
    \fn PDF_Base * PDF_MRST99::GetCopy()
    A method to initialize another electron PDF as exact copy of the current one.
    This is needed for the initial state shower of APACIC.
  */

}  

#endif
