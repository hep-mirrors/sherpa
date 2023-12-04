#ifndef ALARIC__Shower__Kernel_Key_H
#define ALARIC__Shower__Kernel_Key_H

#include <string>

#include "ATOOLS/Phys/Flavour.H"

namespace MODEL {
  class Single_Vertex;
}

namespace ALARIC {

  class Kernel;

  struct Kernel_Key {

    int m_mode, m_swap, m_type;
    std::string m_lfid;
    MODEL::Single_Vertex *p_v;
    ATOOLS::Flavour_Vector m_fl;
    Kernel *p_k;

    Kernel_Key(MODEL::Single_Vertex *const v,
	       const int mode,const int swap,
	       const int type,const std::string &lfid):
      p_k(NULL), m_mode(mode), m_swap(swap), m_type(type),
      m_lfid(lfid), p_v(v)
    {}
    Kernel_Key(const ATOOLS::Flavour_Vector &fl,
	       const int mode,const int type,
	       const std::string &lfid):
      p_k(NULL), m_mode(mode), m_swap(0), m_type(type),
      m_lfid(lfid), p_v(NULL), m_fl(fl)
    {}

  };// end of struct Kernel_Key

}// end of namespace ALARIC

#endif
