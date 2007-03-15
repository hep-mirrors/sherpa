#ifndef Vertex_H
#define Vertex_H

#include "Current.H"

namespace EXTRAXS {

  class Vertex {
  private:

    Current *p_out, *p_a, *p_b;

    bool m_zero;

  public:

    Vertex(Current *const c);

    ~Vertex();

    inline void SetJA(Current *const a) { p_a=a; a->AttachOut(this); }
    inline void SetJB(Current *const b) { p_b=b; b->AttachOut(this); }

    inline void SetZero(const bool zero) { m_zero=zero; }

    inline Current *const JA() const { return p_a; }
    inline Current *const JB() const { return p_b; }

    inline bool Zero() const { return m_zero; }

  };// end of class Vertex

  std::ostream &operator<<(std::ostream &str,const Vertex &v);

}// end of namespace EXTRAXS

#endif
