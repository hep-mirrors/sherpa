#ifndef PDF_LHAPDF_LHAPDF_Fortran_Interface_H
#define PDF_LHAPDF_LHAPDF_Fortran_Interface_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PDF {
  class LHAPDF_Fortran_Interface : public PDF_Base {
  private:
    std::string m_set;
    int         m_member;
    int         m_anti;
    double      m_f[13];
    std::vector<double> m_fv;
  public:
    LHAPDF_Fortran_Interface(const ATOOLS::Flavour,std::string,int);
    PDF_Base * GetCopy();

    void   Calculate(double,double);
    double AlphaSPDF(double);
    double GetXPDF(const ATOOLS::Flavour);

  };


  /*!
    \class LHAPDF_Fortran_Interface
    \brief An interface to the package of various sets of parton distribution functions organized 
    according to the Les Houches accord. (latest update for version LHAPDFv3)

    This class provides an interface to a package of various parton distribution functions
    that are accessible according to the standards put down in the Les Houches accord. The
    underlying idea of this package is to start with sets of pdfs at some scale and to do
    a numerical evolution to a grid of scales. It therefore has to initialize every set that
    is called for the first time which might take some time. Further calls to such a set
    are not so time-consuming, then the grid is already evolved and stored and has just to be
    read-in.

    The package is in Fortran and is provided by Walter Giele, and it allows for an analysis of
    pdf uncertainties according to the ideas of Walter Giele, Stephane A. Keller, David A. Kosower. 
    For more information on their approach, see
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?eprint=hepph-0104052">this publication</A>
    and the <A HREF="http://vircol.fnal.gov/theory.htm">theory section of the LHAPDF manual</A>.
    The package includes the sets:
    - <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?eprint=hepph-0201195">cteq6</A>
    - <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?eprint=hepph-0110215">MRST2001</A>
    - <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?eprint=hepph-0011002">Alekhin</A>
    - <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?eprint=hepph-9912439">Botje</A>

    The Fortran package is also downloadable from <A HREF="http://vircol.fnal.gov/">here</A>.
  */
  /*!
    \var std::string LHAPDF_Fortran_Interface::m_path
    The path from where the initial pdfs and grids are to be loaded in. This will be handed over to
    the corresponding Fortran subroutine as a char*.
  */
  /*!
    \var std::string LHAPDF_Fortran_Interface::m_set
    The name of the set to be used. Available are : cteq6, MRST2001, Alekhin10001, Alekhin1000, 
    Botje100, Botje1000.
  */
  /*!
    \var int LHAPDF_Fortran_Interface::m_member
    The subset of the set named above. Usually "1" is the optimized version with the best fit to data.
  */
  /*!
    \var int LHAPDF_Fortran_Interface::m_anti
    To denote particles and antiparticles, m_anti = 1, -1, respectively.
  */
  /*!
    \var double LHAPDF_Fortran_Interface::m_f[13]
    The array of the 13 pdfs.
  */
  /*!
    \fn LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour,
                                                           std::string,int,std::string,bool &)
    In the constructor the corresponding PDF of the LHAPDF package is initialized. This results
    in either a read-in of an already prepared grid or its calculation.
    The list of partons is filled and the boolean is set to true. this will in the following
    prevent the Fortran code to be initialised twice.
  */
  /*!
    \fn PDF_Base * LHAPDF_Fortran_Interface::GetCopy()
    A method to initialize another electron PDF as exact copy of the current one.
    This is needed for the initial state shower of APACIC.
  */
  /*!
    \var double LHAPDF_Fortran_Interface::AlphaSPDF(double scale2)
    Gives the value of \f$\alpha_s\f$ at a given scale squared that is used internally.
    This might become important to enforce self-consistency of calculations.
  */
  /*!
    void LHAPDF_Fortran_Interface::Calculate(const double, const double)
    Calculates the initialized set at a given energy fraction and a scale. Note that within the SHERPA
    package usually scales are handled as squares of scales, whereas within the LHAPDF package 
    the scales come without squaring.
  */
  /*!
    \fn double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &)
    This returns x times one specific pdf out of the set of all calculated ones - it checks 
    whether the flavour passed over as argument is included in the list of legitimate constituents, 
    PDF_Base::m_partons and selects the appropriate one.
  */
}
 
#endif
