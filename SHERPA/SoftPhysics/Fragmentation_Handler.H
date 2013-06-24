#ifndef SHERPA_SoftPhysics_Fragmentation_Handler_H
#define SHERPA_SoftPhysics_Fragmentation_Handler_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Exception.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif

#ifdef USING__Ahadic
namespace AHADIC {
  class Ahadic;
}
#endif

namespace SHERPA {

  class Fragmentation_Handler: public ATOOLS::Terminator_Object {
  private:

    std::string m_dir,m_file, m_sfile;
    std::string m_fragmentationmodel;
    int         m_mode;
    bool        m_shrink, m_flagpartonics;

#ifdef USING__Ahadic
    AHADIC::Ahadic * p_ahadic;
#endif
#ifdef USING__PYTHIA
    Lund_Interface * p_lund;
#endif

    ATOOLS::Return_Value::code ExtractSinglets(ATOOLS::Blob_List *);
    void Shrink(ATOOLS::Blob_List *);
    void PrepareTerminate(); 
  public:
    Fragmentation_Handler(std::string,std::string);
    ~Fragmentation_Handler();

    ATOOLS::Return_Value::code PerformFragmentation(ATOOLS::Blob_List *,
						    ATOOLS::Particle_List * = NULL);
    
    const std::string &           FragmentationModel() const { return m_fragmentationmodel; }
#ifdef USING__PYTHIA
    inline Lund_Interface * GetLundInterface()   const { return p_lund; }
#endif
#ifdef USING__Ahadic
    inline AHADIC::Ahadic * GetAhadic()          const { return p_ahadic; }
#endif
  };// end of class Fragmentation_Handler

  /*!
    \class Fragmentation_Handler
    \brief The class steering the selection of fragmentation models and their actual work.

    This class is responsible for the steering of fragmentation, i.e. the transition of partons
    into (primary) hadrons. Sooner or later, some models for that will be available, this
    handler will then select and initialize one of them according to the input from the 
    corresponding data file. During event generation, the method PerformFragmentation
    will hand over the task of the actual fragmentation to this particular model. So far, however,
    only the Lund string fragmentation is implemented through an interface to Pythia 6.214,
    the Lund_Fortran_Interface.
  */
  /*!
    \var std::string Fragmentation_Handler::m_dir
    The directory where the file with data to select and steer the actual fragmentation model
    can be found. The default is "./".
  */ 
  /*!
    \var std::string Fragmentation_Handler::m_file
    The file containing the data to select and steer the actual fragmentation model.
    The default is "Fragmentation.dat".
  */ 
  /*!
    \var std::string Fragmentation_Handler::m_fragmentationmodel
    A string naming the fragmentation model to be used. Default (and only implemented
    model so far) is "Lund", i.e. the Lund string model of Pythia, made available through
    the Lund_Fortran_Interface.
  */ 
  /*!
    \var double Fragmentation_Handler::m_lund_a
    This is the parameter \f$a = a_1 = a_2\f$ of the Lund symmetric fragmentation function
    \f[
    f(z)\propto \frac{1}{z} z^{a_1} \left(\frac{1-z}{z}\right)^{a_2}
                \exp\left(-\frac{b m_\perp^2}{z}\right)\,.
    \f]
    This function governs the break-up of the string, \f$z\f$ and \f$m_\perp\f$ are
    parameters describing the kinematics of the quark-anti quark pair produced in
    this break-up. This parameter will finally set PARJ(41) of Pythia.
    The default is set to 0.4 GeV\f$^{-2}\f$.
  */ 
  /*!
    \var double Fragmentation_Handler::m_lund_b
    This is the parameter \f$b\f$ of the Lund symmetric fragmentation function
    \f[
    f(z)\propto \frac{1}{z} z^{a_1} \left(\frac{1-z}{z}\right)^{a_2}
                \exp\left(-\frac{b m_\perp^2}{z}\right)\,.
    \f]
    This function governs the break-up of the string, \f$z\f$ and \f$m_\perp\f$ are
    parameters describing the kinematics of the quark-anti quark pair produced in
    this break-up. This parameter will finally set PARJ(41) of Pythia.
    The default is set to 0.85 GeV\f$^{-2}\f$.
  */ 
  /*!
    \var double Fragmentation_Handler::m_lund_sigma
    The Gaussian width that governs the transverse distribution of primary hadrons.
    The default is set to 0.36 GeV.
  */ 
  /*!
    \var int Fragmentation_Handler::m_mode
    The actual mode labelling the fragmentation model to be used. At the moment, only
    1, i.e. Lund string fragmentation, and 0, i.e. off is available.
  */ 
  /*! 
    \var Lund_Fortran_Interface * Fragmentation_Handler::p_lund 
    Interface to the Lund String. More checks are needed for hadronic collisions.
  */
  /*!
    \fn bool Fragmentation_Handler::ExtractSinglets(ATOOLS::Blob_List *,ATOOLS::Parton_List *)
    This method is used to extract colour ordered chains of partons before hadronization.
    It loops over all blobs in the list and over all outgoing partons in these blobs, until
    an active parton in the colour triplet state (i.e. first colour index different, second
    colour index equal to zero) has been found. This parton is fed into FindConnected
    which will construct recursively the full chain spanned by this parton. All partons
    absorbed in these chains are set to decayed, i.e. their status flag is switched to 2.
    ExtractSinglets will continue finding active triplets and constructing chains until all
    active blobs of the list are spent.
   */
  /*!
    \fn bool Fragmentation_Handler::FindConnected(ATOOLS::Blob_List *,ATOOLS::Parton *,
			                          ATOOLS::Blob *)
    Used for the recursive construction of singlet chains. This method compares the first
    (triplet) colour index of the original parton with the second (anti-triplet) colour index
    of all suitable (i.e. outgoing) active partons of the blob list until a match is found. 
    Its status is set to 2, "decayed". If the match ends the colour chain, in other words, 
    if it is in a colour anti-triplet state (first component equals zero) true is returned. If 
    the matching parton has a triplet index (in other words, if it is in a colour octet state), 
    it serves as original parton for the recursive call to FindConnected.
  */
  /*!
    \fn Fragmentation_Handler::Fragmentation_Handler(std::string,std::string)
    In the constructor the fragmentation model is chosen and its parameters are read in,
    before the according specific model is initialized. Currently, only the Lund string
    is accessible.
  */
  /*!
    \fn Fragmentation_Handler::~Fragmentation_Handler()
    Deletes the fragmentation models.
  */
  /*!
    \fn bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List *,
                                                         ATOOLS::Parton_List * = NULL)
    After extracting singlet chains from the blob list and filling them as incoming partons 
    into blobs of the type "Fragmentation" through the method  ExtractSinglets, these
    new blobs are boosted into their c.m. frame, treated through the Hadronize method of the 
    according fragmentation model, producing hadrons, and then boosted back into the lab 
    system.

    \todo I'm not sure about the boosts here. What happens to the hadron decay blobs ?
  */
  /*!
    \fn Lund_Fortran_Interface * Fragmentation_Handler::GetLundFortranInterface()
    Returns the pointer to te Lund_Fortran_Interface, p_lund.
  */
  /*!
    \fn std::string Fragmentation_Handler::FragmentationModel()
    Returns the label of the fragmentation model, m_fragmentationmodel.
  */
}
#endif

