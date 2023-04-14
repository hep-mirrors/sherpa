#ifndef SHERPA_Soft_Photons_H
#define SHERPA_Soft_Photons_H

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
  class Scoped_Settings;
}

// namespace PHASIC {
//   struct Process_Info;
//   struct Subprocess_Info;
// }

// namespace MODEL {
//   class Single_Vertex;
//   typedef std::vector<Single_Vertex *> Vertex_List;
// }

namespace SHERPA {

  class Soft_Photons: public Event_Phase_Handler {
  private :
    bool                      m_on;
    Matrix_Element_Handler  * p_mehandler;
    YFS::YFS_Handler        * p_yfshandler;
    ATOOLS::NLO_subevtlist  * p_newsublist;

    // void RegisterDefaults(ATOOLS::Scoped_Settings&);

  public :

    // constructor
    Soft_Photons(Matrix_Element_Handler*);

    // destructor
    ~Soft_Photons();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);

  };// end of class Soft_Photons

}
#endif

