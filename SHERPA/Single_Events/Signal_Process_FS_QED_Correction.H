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
  typedef std::vector<Particle* > Particle_Vector;
}

namespace PHASIC {
  struct Process_Info;
  struct Subprocess_Info;
}

namespace SHERPA {

  typedef std::vector<ATOOLS::Blob *> Blob_Vector;
  typedef std::vector<const PHASIC::Subprocess_Info *> SubInfoVector;

  class Signal_Process_FS_QED_Correction: public Event_Phase_Handler {
  private :
    bool                      m_on;
    bool                      m_qed;
    Matrix_Element_Handler  * p_mehandler;
    Soft_Photon_Handler     * p_sphotons;

    std::map<std::string,SubInfoVector> m_proc_lep_map;

    bool            PutOnMassShell(const ATOOLS::Particle_Vector&);
    ATOOLS::Flavour DetermineResonanceFlavour(const ATOOLS::Particle_Vector&);
    ATOOLS::Vec4D   MomentumSum(const ATOOLS::Particle_Vector&);

    void FindSubProcessInfosContainingLeptons(const PHASIC::Process_Info&,
                                              SubInfoVector&);
    void FindSubProcessInfosContainingLeptons(const PHASIC::Subprocess_Info&,
                                              SubInfoVector&);

    Blob_Vector BuildResonantBlobs(ATOOLS::Particle_Vector&);

    bool    ContainsNoAmbiguities(const ATOOLS::Particle_Vector&);
    void    FillBlob(ATOOLS::Blob *,
                     const PHASIC::Subprocess_Info&,
                     ATOOLS::Particle_Vector&);
    void    FillBlob(ATOOLS::Blob *,
                     const ATOOLS::Flavour&,
                     ATOOLS::Particle_Vector&);
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

