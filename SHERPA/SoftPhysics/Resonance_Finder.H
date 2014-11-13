#ifndef SHERPA_SoftPhysics_Resonance_Finder_H
#define SHERPA_SoftPhysics_Resonance_Finder_H

#include <map>
#include <string>
#include <vector>
#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {
  class Blob;
  class Data_Reader;
  class Flavour;
  class Particle;
  typedef std::vector<Blob *> Blob_Vector;
  typedef std::vector<Flavour> Flavour_Vector;
  typedef std::vector<Particle *> Particle_Vector;
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
  typedef std::vector<const PHASIC::Subprocess_Info *> SubInfoVector;
  class Matrix_Element_Handler;

  class Resonance_Finder {
  private:
    bool   m_on;
    double m_resdist;
    Matrix_Element_Handler  * p_mehandler;
    std::map<std::string,SubInfoVector>      m_proc_lep_map;
    std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;

    void ScanModelForEWResonances();
    void FindProcessPossibleResonances(const ATOOLS::Flavour_Vector& fv,
                                       MODEL::Vertex_List& vlist);
    void InitialiseHelperParticles();
    void IdentifyEWSubprocesses();

    void FindSubProcessInfosContainingLeptons(const PHASIC::Process_Info&,
                                              SubInfoVector&);
    void FindSubProcessInfosContainingLeptons(const PHASIC::Subprocess_Info&,
                                              SubInfoVector&);

    void FillBlob(ATOOLS::Blob *,
                  const PHASIC::Subprocess_Info&,
                  ATOOLS::Particle_Vector&);
    void FillBlob(ATOOLS::Blob *,
                  const ATOOLS::Flavour&,
                  ATOOLS::Particle_Vector&);

    bool FindResonances(ATOOLS::Particle_Vector&,
                        std::vector<ATOOLS::Particle_Vector>&,
                        ATOOLS::Flavour_Vector&,
                        const MODEL::Vertex_List&);

    bool ContainsNoAmbiguities(const ATOOLS::Particle_Vector&);

    ATOOLS::Flavour DetermineGenericResonance(const ATOOLS::Particle_Vector&);
    ATOOLS::Vec4D   MomentumSum(const ATOOLS::Particle_Vector&);

    public:
    Resonance_Finder(ATOOLS::Data_Reader* reader,
                     Matrix_Element_Handler * meh);
    Resonance_Finder();
    ~Resonance_Finder();

    void BuildResonantBlobs(ATOOLS::Particle_Vector&,
                            ATOOLS::Blob_Vector& blobs);

  };
}// end of namespace SHERPA

#endif
