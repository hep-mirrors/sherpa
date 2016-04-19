#ifndef METOOLS_Explicit_Vertex_Key_H
#define METOOLS_Explicit_Vertex_Key_H

#include "METOOLS/Explicit/Dipole_Info.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/STL_Tools.H"

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
  private:

    Vertex_Key
    (const std::vector<Current*> &j,
     Current *const c,MODEL::Model_Base *const model,
     MODEL::Single_Vertex *const mv=NULL,const std::string &p="",
     Vertex *const v=NULL,Color_Calculator *const cc=NULL,
     Lorentz_Calculator *const lc=NULL);

  public:

    Current *p_c, *p_k, *p_kt;
    std::vector<Current*> m_j;

    MODEL::Model_Base    *p_model;
    MODEL::Single_Vertex *p_mv;

    std::string m_p;
    size_t      m_n;
    mutable int m_d;

    mutable std::string m_id;

    mutable Vertex *p_v;

    mutable Color_Calculator   *p_cc;
    mutable Lorentz_Calculator *p_lc;

    Dipole_Info *p_dinfo;

    static ATOOLS::AutoDelete_Vector<Vertex_Key> s_objects;

    static Vertex_Key *New
    (const std::vector<Current*> &j,
     Current *const c,MODEL::Model_Base *const model,
     MODEL::Single_Vertex *const mv=NULL,const std::string &p="",
     Vertex *const v=NULL,Color_Calculator *const cc=NULL,
     Lorentz_Calculator *const lc=NULL);

    void Delete();

    std::string Type() const;

    const std::string &ID() const;

    ATOOLS::Flavour Fl(const size_t &i) const;

    Current *J(const size_t &i) const;

  };// end of struct Vertex_Key

}// end of namespace METOOLS

#endif
