#ifndef AHADIC_Tools_Hadron_Init_H
#define AHADIC_Tools_Hadron_Init_H

#include "ATOOLS/Org/Data_Reader.H"

namespace AHADIC {

 class Hadron_Init {
  public:
    Hadron_Init() {}
    void Init();
    void OverWriteProperties(ATOOLS::Data_Reader& overwrite);
  };

}

#endif
