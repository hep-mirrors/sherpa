#ifndef METOOLS_Explicit_Vertex_Key_H
#define METOOLS_Explicit_Vertex_Key_H

#include "METOOLS/Explicit/Dipole_Info.H"
#include "ATOOLS/Phys/Flavour.H"

#include <string>

namespace MODEL {
  class Model_Base;
  class Single_Vertex;
}

namespace METOOLS {

  class Current;
  class Vertex;
  class Color_Calculator;
  class Lorentz_Calculator;

  struct Vertex_Key {

    Current  *p_a, *p_b, *p_c, *p_e, *p_k, *p_kt;

    MODEL::Model_Base    *p_model;
    MODEL::Single_Vertex *p_mv;

    std::string m_p;
    size_t      m_n;
    mutable int m_d;

    mutable Vertex *p_v;

    mutable Color_Calculator   *p_cc;
    mutable Lorentz_Calculator *p_lc;

    Dipole_Info *p_dinfo;

    Vertex_Key
    (Current *const a,Current *const b,Current *const e,
     Current *const c,MODEL::Model_Base *const model,
     MODEL::Single_Vertex *const mv=NULL,const std::string &p="",
     Vertex *const v=NULL,Color_Calculator *const cc=NULL,
     Lorentz_Calculator *const lc=NULL);

    std::string Type() const;

    std::string ID(const int dir=0) const;

    ATOOLS::Flavour FlA() const;
    ATOOLS::Flavour FlB() const;

    Vertex_Key SwapAB() const;
    Vertex_Key SwapBE() const;
    Vertex_Key SwapEA() const;

    Current *J(const size_t &i) const;

    bool operator<(const Vertex_Key &k) const;

  };// end of struct Vertex_Key

}// end of namespace METOOLS

#endif
