#ifndef Hadron_Init_H
#define Hadron_Init_H

#include "Data_Reader.H"

namespace AHADIC {

 class Hadron_Init {
  public:
    Hadron_Init() {}
    void Init();
    void OverWriteProperties(ATOOLS::Data_Reader& overwrite);
  };

}

#endif
