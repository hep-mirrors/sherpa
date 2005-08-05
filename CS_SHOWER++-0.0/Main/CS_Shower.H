#ifndef CS_Shower_h
#define CS_Shower_h

#include "Shower.H"
#include "Singlet.H"
#include "Blob_List.H"
#include "Blob.H"

namespace CS_SHOWER {
  class CS_Shower {
  private:
    Shower         m_shower;
    All_Singlets * p_allsinglets;

    void ExtractPartons(ATOOLS::Blob *);
  public:
    CS_Shower();
    ~CS_Shower();

    bool PerformShowers(All_Singlets * =NULL);
    bool ExtractPartons(ATOOLS::Blob_List *);
    void PrepareAllSinglets();

    inline All_Singlets * GetAllSinglets() { return p_allsinglets; }
  };
}

#endif
