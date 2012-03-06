#ifndef SHERPA_Single_Events_Jet_Evolution_H
#define SHERPA_Single_Events_Jet_Evolution_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"

namespace SHERPA {
  
  class Hard_Decay_Handler;

  class Jet_Evolution: public Event_Phase_Handler {
  private :
    PertInterfaceMap          m_interfaces;

    ATOOLS::Return_Value::code AttachShowers(ATOOLS::Blob *,ATOOLS::Blob_List *,
					     Perturbative_Interface *);
    
    void Reset();
    void Update(const ATOOLS::Blob *,const size_t,
                Perturbative_Interface *interface);
    bool DefineInitialConditions(const ATOOLS::Blob *,const ATOOLS::Blob_List *,
                                 Perturbative_Interface *interface);
    void AftermathOfNoShower(ATOOLS::Blob *,ATOOLS::Blob_List *);
    void AftermathOfSuccessfulShower(ATOOLS::Blob *,ATOOLS::Blob_List *,
				     Perturbative_Interface *);    
  public :

    // constructor
    Jet_Evolution(Matrix_Element_Handler *, Hard_Decay_Handler*,
                  Decay_Handler_Base *,MI_Handler *, Soft_Collision_Handler *,
                  const Shower_Handler_Map &);

    // destructor
    ~Jet_Evolution();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);

  };// end of class Jet_Evolution

  /*!
    \file 
    \brief Contains the class SHERPA::Jet_Evolution
  */
  /*!
    \class Jet_Evolution
    \brief The handler for the perturbative shower evolution around hard production processes or decays.

    This class deals with the perturbative evolution of a parton ensemble after a hard production process
    or decay has taken place. This evolution is modelled through multiple readiation of secondaries,
    formulated in a shower language. This class therefore contains a pointer to a Shower_Handler, to
    accomodate a suitable merging of the regimes of jet production (as modelled, e.g., by the Signal_Processes)
    and jet evolution, i.e. the showers, an interface has to be defined. This interface will ensure
    - a proper filling of the phase space for the emission of secondaries without double counting (over or 
    under counting) and it depends to a large extend on the actual shower and the information passed
    by the Matrix_Element_Handler. Consequently, this handler also has to access such information which is
    achieved by a pointer to the actual Matrix_Element_Handler.
  */
  /*!
    \var Matrix_Element_Handler * Jet_Evolution::p_mehandler
    The pointer to the Matrix_Element_Handler of choice. This choice is also reflected in the actual
    realization of the matrix element - parton shower interface (the Perturbative_Interface). At the moment,
    only two Matrix_Element_Handler's are available, namely Amegic and the internal 2->2 matrix elements which
    are in fact quite limited in scope.
  */
  /*!
    \var Shower_Handler * Jet_Evolution::p_showerhandler
    The pointer to the Shower_Handler of choice. rightnow, only a Pythia-type parton shower (i.e. ordering by 
    virtual mass plus angular ordering imposed "by hand" to ensure coherence) is available, in the near
    future, however, there will be also a dipole shower available.
  */
  /*!
    \var Perturbative_Interface * Jet_Evolution::p_interface
    The pointer to the interface between matrix element and parton/dipole shower. This interface basically
    is responsible for taking the blob stemming from the hard production process or decay and translate the
    corresponding partons into a parton history suitable for showering (in the actual case of a parton shower
    this boils down to the fillnig of the partons into a binary tree structure of 1->2 decays of "mothers"
    into "daughters"). It should be stressed that every combination of matrix element generator and shower
    needs a specific realization of this due to different information obtainalble from the matrix elements and
    different shower variables used in the shower.
  */
  /*!
    \fn Jet_Evolution::Jet_Evolution(Matrix_Element_Handler *,Shower_Handler *)
    In the constructor the pointers to the Matrix_Element_Handler and the Shower_Handler are copied. Also, 
    the constructor sets the name and type of the Event_Phase_Handler, m_name and m_type. The name is 
    "Jet_Evolution : " plus the name of the shower generator, and the type is "Perturbative". Then a suitable
    interface between the matrix elements and the parton shower is instantiated. At the moment, only two
    such interfaces are available, namely Amegic_Apacic_Interface and SimpleXS_Apacic_Interface.
  */
  /*!
    \fn Jet_Evolution::~Jet_Evolution()
    Deletes the Perturbative_Interface.
  */
  /*!
    \fn bool Jet_Evolution::Treat(ATOOLS::Blob_List *, double &)
    In this method all blobs with status equal to one and of type "Signal Process :" are dealt with by
    the follwonig steps:
    \todo Andreas, we have to fill this.
  */
  /*!
    \fn void Jet_Evolution::CleanUp()
    Cleans up the shower handler by its corresponding method.
  */
}
#endif

