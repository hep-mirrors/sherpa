#ifndef CTEQ6_Fortran_Interface_H
#define CTEQ6_Fortran_Interface_H

#include "PDF_Base.H"
#include "Flavour.H"

namespace PDF {
  class CTEQ6_Fortran_Interface : public PDF_Base {
  private:
    std::string m_set,m_path;
    int         m_member;
    int         m_anti;
    double      m_f[11];
  public:
    CTEQ6_Fortran_Interface(const APHYTOOLS::Flavour,std::string,int,std::string);
    PDF_Base * GetCopy();
    void   Output();

    void   Calculate(const double, const double);
    double AlphaSPDF(double);
    double GetXPDF(const APHYTOOLS::Flavour &);

    double GetXMin()   { return 0.; }
    double GetXMax()   { return 1.; }
    double GetQ2Min()  { return 1.; }
    double GetQ2Max()  { return 1.e12; }
    APHYTOOLS::Flavour Bunch()                { return m_bunch; }
    std::vector<APHYTOOLS::Flavour> Partons() { return m_partons; }
  };


  /*!
    \class CTEQ6_Fortran_Interface
    \brief An interface to the cteq6 PDF package.
  */
  /*!
    \var std::string CTEQ6_Fortran_Interface::m_path
    The path from where the initial pdfs and grids are to be loaded in. This will be handed over to
    the corresponding Fortran subroutine as a char*.
  */
  /*!
    \var std::string CTEQ6_Fortran_Interface::m_set
    The name of the set to be used. Available are : 

    1    CTEQ6M   Standard MSbar scheme   cteq6m.tbl
    2    CTEQ6D   Standard DIS scheme     cteq6d.tbl
    3    CTEQ6L   Leading Order           cteq6l.tbl
    4    CTEQ6L1  Leading Order           cteq6l1.tbl

  */
  /*!
    \var int CTEQ6_Fortran_Interface::m_member
  */
  /*!
    \var double CTEQ6_Fortran_Interface::m_f[13]
    The array of the 13 pdfs.
  */
  /*!
    \var double CTEQ6_Fortran_Interface::AlphaSPDF(double scale2)
    Gives the value of \f$\alpha_s\f$ at a given scale squared that is used internally.
    This might become important to enforce self-consistency of calculations.
  */
  /*!
    void CTEQ6_Fortran_Interface::Calculate(const double, const double)
    Calculates the initialized set at a given energy fraction and a scale. Note that within the SHERPA
    package usually scales are handled as squares of scales, whereas within the CTEQ6 package 
    the scales come without squaring.
  */
}
 
#endif
