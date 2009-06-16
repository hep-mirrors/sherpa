#ifndef HADRONS_Main_Mixing_Handler_H
#define HADRONS_Main_Mixing_Handler_H

namespace ATOOLS {
  class Blob;
  class Blob_List;
  class Particle;
}

namespace HADRONS {
  class Hadron_Decay_Map;
  class Hadron_Decay_Table;

  class Mixing_Handler {
    private:
      Hadron_Decay_Map* p_decaymap;
      double DetermineMixingTime(ATOOLS::Particle* decayer, bool checkforpartstatus);

    public:
      Mixing_Handler(Hadron_Decay_Map* decaymap);
      ~Mixing_Handler();
      bool PerformMixing(ATOOLS::Particle* decayer, ATOOLS::Blob_List* bloblist);
      bool SetCPAsymmetries(ATOOLS::Particle* decayer, Hadron_Decay_Table* table);
      void ResetCPAsymmetries(Hadron_Decay_Table* table);
  };
}

#endif
