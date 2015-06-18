#ifndef DIRE__Shower__Kernel_Key_H
#define DIRE__Shower__Kernel_Key_H

#include <string>

namespace ATOOLS {
  class Data_Reader;
}

namespace MODEL {
  class Single_Vertex;
}

namespace DIRE {

  class Kernel;

  struct Kernel_Key {

    int m_mode, m_type;
    MODEL::Single_Vertex *p_v;
    Kernel *p_k;
    ATOOLS::Data_Reader *p_rd;

    Kernel_Key(MODEL::Single_Vertex *const v,
	       const int mode,const int type,
	       ATOOLS::Data_Reader *const rd):
      p_k(NULL), m_mode(mode), m_type(type), p_v(v),
      p_rd(rd) {}

  };// end of struct Kernel_Key

}// end of namespace DIRE

#endif
