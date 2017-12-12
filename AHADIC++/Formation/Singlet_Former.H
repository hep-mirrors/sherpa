#ifndef AHADIC_Formations_Singlet_Former_H
#define AHADIC_Formations_Singlet_Former_H

#include "AHADIC++/Tools/Singlet_Tools.H"
#include "ATOOLS/Phys/Blob.H"
#include <list>

namespace AHADIC {
  class Singlet_Former {
  private:
    std::list<Singlet *>        * p_singlets;
    std::list<ATOOLS::Particle *> m_colparts;
    
    void ExtractOutgoingCols(ATOOLS::Blob * blob);
    Singlet          * MakeAnother();
    ATOOLS::Particle * FindStart();
    void FormSinglets();
  public:
    Singlet_Former(std::list<Singlet *> * singlets);
    ~Singlet_Former();
    bool Extract(ATOOLS::Blob * blob);
  };
}

#endif
