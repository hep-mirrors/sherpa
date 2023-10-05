#ifndef SHERPA_SoftPhysics_Singlet_Sorter_H
#define SHERPA_SoftPhysics_Singlet_Sorter_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include <vector>

namespace SHERPA {
  class Singlet_Sorter {
  private :
    std::list<ATOOLS::Particle *>   m_hadrons;
    ATOOLS::Part_List             * p_partlist;
    std::list<ATOOLS::Part_List * > m_partlists;
    ATOOLS:: Momenta_Stretcher      m_stretcher;
    long unsigned int               m_calls, m_fails;

    bool HarvestParticles(ATOOLS::Blob_List *);
    bool FillParticleLists(ATOOLS::Blob * blob);
    void DealWithHadrons(ATOOLS::Blob_List * bloblist);
    bool DecomposeIntoSinglets();
    bool NextSinglet(ATOOLS::Part_List & sorted,const bool trip);
    ATOOLS::Particle  * FindNext(const size_t col);
  public:
    Singlet_Sorter();
    ~Singlet_Sorter();
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *);
    ATOOLS::Blob * MakeBlob();
    void ResetPartLists();
  };
}

#endif
