#ifndef Hdecay_Fortran_Interface_h
#define Hdecay_Fortran_Interface_h

#include "Spectrum_Generator_Base.H"
#include "Data_Read.H"

namespace HDECAY {
  class Hdecay_Fortran_Interface : public MODEL::Spectrum_Generator_Base {
  private :
    double hmass;
    double brb,brtau,brmu,brc,brs,brt;
    double brw,brz,brgg,brgamma,brzgamma,hwidth;

    void   CalculateEffectiveCouplings(std::string);
  protected :

  public :
    Hdecay_Fortran_Interface(AORGTOOLS::Data_Read *,MODEL::Model_Base *);
    ~Hdecay_Fortran_Interface() {};
    void Run(std::string);
    void FillMasses() {}
    void FillDecays() {}
  };

  /*!
    \class Hdecay_Fortran_Interface
    \brief A wrapper around Hdecay.

    This is a wrapper around the Fortran code
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?rawcmd=FIND+eprint+hep-ph%2F9704448">
    Hdecay</A> by A. Djouadi, J. Kalinowski, and M. Spira. It calculates the branching ratios
    of Higgs bosons in the Standard Model and its minimal supersyhmmetreic extension, the MSSM.

    At the moment, only the SM mode is used and only to calculate the total branching width. In
    the future an option is to be implemented that uses the individual branching ratios to
    define the branching ratios in the program and to define the tree-level couplings of the
    Higgs boson(s) under consideration.

    A warning is in order here : Hdecay performs only, if all fermion masses of the second
    and third generation are larger than 0!
  */
  /*!
    \var double Hdecay_Fortran_Interface::hmass
    The mass of the SM Higgs boson.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brb
    The branching ratio for \f$H\to b\bar b\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brtau
    The branching ratio for \f$H\to \tau\bar\tau\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brmu
    The branching ratio for \f$H\to \mu\bar\mu\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brc
    The branching ratio for \f$H\to c\bar c\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brs
    The branching ratio for \f$H\to s\bar s\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brt
    The branching ratio for \f$H\to t\bar t\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brw
    The branching ratio for \f$H\to W^+W^-\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brz
    The branching ratio for \f$H\to ZZ\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brgg
    The branching ratio for \f$H\to gg\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brgamma
    The branching ratio for \f$H\to \gamma\gamma\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::brzgamma
    The branching ratio for \f$H\to Z\gamma\f$.
  */
  /*!
    \var double Hdecay_Fortran_Interface::hwidth
    The total width of the Higgs boson.
  */
  /*!
    \fn void Hdecay_Fortran_Interface::CalculateEffectiveCouplings(std::string);
    At the moment this method extracts the various branching ratios and the overall width
    of the SM Higgs boson and fills the corresponding parameters and decay channels. In the
    future these decay channels are to have an effect on the actually used values and
    on the tree-level couplings in AMEGIC++.
    \todo Mirror the decay information to the outside and update the vertex couplings.
  */
  /*!
    \fn Hdecay_Fortran_Interface::Hdecay_Fortran_Interface(AORGTOOLS::Data_Read *,MODEL::Model_Base *);
    Copies the pointers to the model that owns this spectrum generator and to the
    AORGTOOLS::Data_Read object.
  */
  /*!
    \fn void Hdecay_Fortran_Interface::Run(std::string);
    Before running the Fortran code, some more parameters are to be read in. in particular these are
    - the effective strange mass, default is \f$m_s = 0.48\,.\f$
    - the number of light flavour in the \f$gg\f$ decay mode, default is 5.
    - whether only on-shell W's and Z's are allowed, default is set to 1. (\todo effect?)
    - two more flags \todo find this out !
  */
}

#endif


