#ifndef DIRE__Shower__Kernel_Key_H
#define DIRE__Shower__Kernel_Key_H

#include <string>

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class Data_Reader;
}

namespace MODEL {
  class Single_Vertex;
}

namespace DIRE {

  class Kernel;

  struct Kernel_Key {

    int m_mode, m_swap, m_type;
    std::string m_lfid;
    MODEL::Single_Vertex *p_v;
    ATOOLS::Flavour_Vector m_fl;
    Kernel *p_k;
    ATOOLS::Data_Reader *p_rd;

    Kernel_Key(MODEL::Single_Vertex *const v,
	       const int mode,const int swap,const int type,
	       ATOOLS::Data_Reader *const rd):
      p_k(NULL), m_mode(mode), m_swap(swap), m_type(type), p_v(v),
      p_rd(rd) {}
    Kernel_Key(const ATOOLS::Flavour_Vector &fl,
	       const int mode,const int type,
	       ATOOLS::Data_Reader *const rd,
	       const std::string &lfid):
      p_k(NULL), m_mode(mode), m_type(type),
      m_lfid(lfid), p_v(NULL), m_fl(fl), p_rd(rd) {}

  };// end of struct Kernel_Key

}// end of namespace DIRE

#endif
