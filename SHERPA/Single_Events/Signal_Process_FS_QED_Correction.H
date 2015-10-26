#ifndef SHERPA_Single_Events_Signal_Process_FS_QED_Correction_H
#define SHERPA_Single_Events_Signal_Process_FS_QED_Correction_H

#include "ATOOLS/Math/Vector.H"
#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

namespace ATOOLS {
  class Blob;
  class Flavour;
  class Particle;
  typedef std::vector<Blob *> Blob_Vector;
  typedef std::vector<Particle* > Particle_Vector;
}

namespace PHASIC {
  struct Process_Info;
  struct Subprocess_Info;
}

namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *> Vertex_List;
}

namespace SHERPA {

  class Signal_Process_FS_QED_Correction: public Event_Phase_Handler {
  private :
    bool                      m_on;
    bool                      m_qed;
    bool                      m_forceqedonnloqcd;
    Matrix_Element_Handler  * p_mehandler;
    Soft_Photon_Handler     * p_sphotons;
    ATOOLS::NLO_subevtlist  * p_newsublist;

    bool            PutOnMassShell(const ATOOLS::Particle_Vector&);
    bool            ModifySubEvtList(ATOOLS::Blob *,
                                     const ATOOLS::Particle_Vector&,
                                     const ATOOLS::Blob_Vector&);
    void            DeleteNewSubList();

  public :

    // constructor
    Signal_Process_FS_QED_Correction(Matrix_Element_Handler *,Soft_Photon_Handler *);

    // destructor
    ~Signal_Process_FS_QED_Correction();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);

  };// end of class Signal_Process_FS_QED_Correction

}
#endif

