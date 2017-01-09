#ifndef AHADIC_Tools_Hadron_Init_H
#define AHADIC_Tools_Hadron_Init_H

namespace ATOOLS {
  class Default_Reader;
}

namespace AHADIC {

 class Hadron_Init {
  public:
    Hadron_Init() {}
    void Init();
    void OverWriteProperties(ATOOLS::Default_Reader&);
  };

}

#endif
