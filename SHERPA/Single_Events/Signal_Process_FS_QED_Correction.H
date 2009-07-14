#ifndef SHERPA_Single_Events_Signal_Process_FS_QED_Correction_H
#define SHERPA_Single_Events_Signal_Process_FS_QED_Correction_H

#include "ATOOLS/Math/Vector.H"
#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

namespace ATOOLS {
  class Flavour;
}

namespace SHERPA {

  class Signal_Process_FS_QED_Correction: public Event_Phase_Handler {
  private :
    bool                      m_on;
    MEHandlersMap           * p_mehandlers;
    Soft_Photon_Handler     * p_sphotons;

    bool    PutOnMassShell(const Particle_Vector&);
    Flavour DetermineResonanceFlavour(const Particle_Vector&);
    Vec4D   MomentumSum(const Particle_Vector&);
  public :

    // constructor
    Signal_Process_FS_QED_Correction(MEHandlersMap *,Soft_Photon_Handler *);

    // destructor
    ~Signal_Process_FS_QED_Correction();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);

  };// end of class Signal_Process_FS_QED_Correction

}
#endif

