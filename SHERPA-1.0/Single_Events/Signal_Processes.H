#ifndef Signal_Processes_h
#define Signal_Processes_h

#include "Event_Phase_Handler.H"
#include "Matrix_Element_Handler.H"
#include "Hard_Decay_Handler.H"
#include "Remnant_Base.H"

namespace SHERPA {
  class Signal_Processes : public Event_Phase_Handler {
  private :
    Matrix_Element_Handler * p_mehandler;
    Hard_Decay_Handler     * p_hdhandler;
    Remnant_Base           * p_remnants[2];
    void FillBlob(ATOOLS::Blob *, const double, const int);
  public :
    Signal_Processes(Matrix_Element_Handler *,Hard_Decay_Handler * = NULL);
    ~Signal_Processes();
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);
  };

  /*!
    \class Signal_Processes
    \brief This is the Event_Phase_Handler concerned with the generation of the hardest
           sub-process in an event, which is called the signal process.

    This class is concerned with the generation of the hardest sub-process in each event through
    the Matrix_Element_Handler that was specified before in the initialization. The idea here is that 
    a blob list containing exactly one blob is passed into this class. This blob is then dealt with
    accordingly by the methods Treat(ATOOLS::Blob_List *, double &) and FillBlob().

    In principle with this one could have various hard interactions of different pairs of incoming 
    particles as well, but at the moment only the possibility for one such hard subprocess per bunch 
    crossing is implemented. 
    \todo Enable more hard sub-processes per bunch crossing.
  */
  /*!
    \var Matrix_Element_Handler *  Signal_Processes::p_mehandler
    The pointer to the Matrix_Element_Handler of choice. This choice is reflected in the name
    of the Handler, namely "Signal_Processes : " plus the name of the Matrix_Element_Handler. At the moment,
    only two such handlers are available, namely Amegic and the internal 2->2 matrix elements which
    are in fact quite limited in scope.

    \todo In the future there should be implemented some citation system for codes and methods used 
          with Sherpa.
    \todo Expand the internal list of processes.
  */
  /*!
    \fn void Signal_Processes::FillBlob(ATOOLS::Blob *)
    Filling the blob means:
    - Setting the position (to the origin),
    - setting its type (actually, this should be made prettier) and its status (to 1)
    - setting its c.m. momentum as defined by the incoming partons
    - translating the list of flavours and momenta stemming from the Matrix_Element_Handler into
      newly instantiated partons with suitable status, info code, and with pointers to production and 
      decay blob, depending on whether they are incoming or outgoing.
  */
  /*!
    \fn Signal_Processes::Signal_Processes(Matrix_Element_Handler *)
    In the constructor the pointer to the Matrix_Element_Handler is copied. Also, the constructor
    sets the name and type of the Event_Phase_Handler, m_name and m_type. The name is "Signal_Processes : " 
    plus the name of the Matrix_Element_Handler, and the type is "Perturbative".
  */
  /*!
    \fn Signal_Processes::~Signal_Processes()
    Deletes the Matrix_Element_Handler. This might lead to trouble.
  */
  /*!
    \fn bool Signal_Processes::Treat(ATOOLS::Blob_List *, double &)
    The blob list is checked if it contains exactly one blob. This starting blob must have status
    equal to 0 or -1 and it must be of the type "Signal Process : ". Then this Event_Phase_Handler
    will come into action (it should be noted however that structurally more blobs can be treated after
    minor adjustments in this handler and more thinking around it). If the status of the blob
    equals 0, then the Signal_Processes will employ the Matrix_Element_Handler::GenerateOneEvent() -method
    to generate a set of momenta and flavours for a newly picked process. If the status is -1, then
    the Matrix_Element_Handler will jsut porduce another kinematical configuration for the same
    process through the Matrix_Element_Handler::GenerateSameEvent() -method. The latter case applies, 
    for instance, if the jet evolution throws a veto because of migration of an n-jet configuration
    to an (n-1)-jet configuration when the standard merging prescription for MEs and the parton shower
    is applied. In any case, however GenerateOneEvent() or GenerateSameEvent() is employed,
    the blob in question will be filled accordingly through the FillBlob(Blob *) method.
  */
  /*!
    \fn void Signal_Processes::CleanUp()
    Nothing has to be done in this Event_Phase_Handler at the moment. However, it is worth thinking about
    shifting some resetting of the internal blob to be filled into this method.
  */
}
#endif

