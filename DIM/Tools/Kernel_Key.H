#ifndef DIM__Shower__Kernel_Key_H
#define DIM__Shower__Kernel_Key_H

#include <string>

namespace MODEL {
  class Single_Vertex;
}

namespace DIM {

  class Kernel;

  struct Kernel_Key {

    int m_mode, m_type;
    MODEL::Single_Vertex *p_v;
    Kernel *p_k;

    Kernel_Key(MODEL::Single_Vertex *const v,
	       const int mode,const int type):
      p_k(NULL), m_mode(mode), m_type(type), p_v(v)
    {}

  };// end of struct Kernel_Key

}// end of namespace DIM

#endif
