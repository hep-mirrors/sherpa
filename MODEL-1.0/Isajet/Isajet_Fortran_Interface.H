#ifndef Isajet_Fortran_Interface_h
#define Isajet_Fortran_Interface_h

#include "Spectrum_Generator_Base.H"
#include "Data_Read.H"
#include "Matrix.H"

namespace ISAJET { 
  class Isajet_Fortran_Interface : public MODEL::Spectrum_Generator_Base {
    void CharginoMasses();
    void NeutralinoMasses();
    void GluinoMasses();
    void HiggsMasses();
    void sNeutrinoMasses();
    void sUpMasses();
    void sDownMasses();
    void sLeptonMasses();

  public:
    Isajet_Fortran_Interface(ATOOLS::Data_Read *,MODEL::Model_Base *);
    ~Isajet_Fortran_Interface();
    void Run(std::string);
    void FillMasses();
    void FillDecays();
  };

  /*!
    \class Isajet_Fortran_Interface
    \brief A wrapper around Isajet

    This is a wrapper around the version 7.67 of the Fortran code
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?rawcmd=FIND+eprint+hep-ph%2F0001086">
    Isajet</A> by H. Baer, F. E. Paige, S. D. Protopopescu, and X. Tata.

    At the moment, only the mass spectra of SUSY is calculated. In the future an option is to be 
    implemented that uses the individual branching ratios to define the branching ratios in the 
    program and to define the corresponding tree-level couplings.
   */
  /*!
    \fn void Isajet_Fortran_Interface::CharginoMasses()
    Reads in the chargino masses, sets them in the Flavour class, reads in the parameters for
    the two mixing matrices \f$Z^\pm\f$ and adds them to the list of complex matrices in 
    the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::NeutralinoMasses()
    Reads in the neutralino masses, sets them in the Flavour class, reads in the parameters for
    the mixing matrix \f$Z^N\f$ and adds it to the list of complex matrices in the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::GluinoMasses()
    Reads in the gluino mass and sets it.
  */
  /*!
    \fn void Isajet_Fortran_Interface::HiggsMasses()
    Reads in the higgs boson masses, sets them in the Flavour class, reads in the parameters for
    the mixing matrices \f$Z^{H,R}\f$ and adds them to the list of complex matrices in the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::sNeutrinoMasses()
    Reads in the masses of the sneutrinos, sets them in the Flavour class, reads in the 
    parameters for the mixing matrix \f$Z^{\nu}\f$ and adds it to the list of complex matrices in 
    the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::sUpMasses()
    Reads in the masses of the up-type squarks, sets them in the Flavour class, reads in the 
    parameters for the mixing matrix \f$Z^{U}\f$ and adds it to the list of complex matrices in 
    the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::sDownMasses()
    Reads in the masses of the down-type squarks, sets them in the Flavour class, reads in the 
    parameters for the mixing matrix \f$Z^{D}\f$ and adds it to the list of complex matrices in 
    the model.
  */
  /*!
    \fn void Isajet_Fortran_Interface::sLeptonMasses()
    Reads in the masses of the sleptons, sets them in the Flavour class, reads in the 
    parameters for the mixing matrix \f$Z^{L}\f$ and adds it to the list of complex matrices in 
    the model.
  */
  /*!
    Isajet_Fortran_Interface(ATOOLS::Data_Read *,MODEL::Model_Base *)
    Copies the pointers to the model that owns this spectrum generator and to the
    ATOOLS::Data_Read object.
  */
  /*!
    \fn void Isajet_Fortran_Interface::Run(std::string)
    Fills in the parameters of the model into suitable arrays and calls isasusy through an 
    interface.
  */
  /*!
    \fn void Isajet_Fortran_Interface::FillMasses()
    This method basically calls all methods xxxMasses in this class.
  */
  /*!
    \fn void Isajet_Fortran_Interface::FillDecays()
    Still empty.
    \todo Fill it!
  */
}
#endif
