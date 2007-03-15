#ifndef Current_H
#define Current_H

#include "Vector.H"
#include "Flavour.H"

#include <vector>

namespace EXTRAXS {

  typedef std::vector<int> Int_Vector;

  class Vertex;

  typedef std::vector<Vertex*> Vertex_Vector;

  struct ct {
    enum type {
      spinor = 1,
      vector = 2,
      tensor = 4
    };
  };// end of struct ct

  std::ostream &operator<<(std::ostream &str,const ct::type &type);

  class CDBG_Amplitude;

  class Current {
  protected:

    CDBG_Amplitude *p_ampl;

    ATOOLS::Flavour m_fl;

    Vertex_Vector m_in, m_out;
    Int_Vector    m_id;

    ATOOLS::Vec4D m_p;

    size_t   m_key;
    ct::type m_type;

  public:

    Current(CDBG_Amplitude *const ampl,const ATOOLS::Flavour &fl);

    virtual ~Current();

    virtual Vertex *GetVertex(Current *const ja,Current *const jb) = 0;

    virtual void SetCurrent() = 0;
    virtual void Evaluate() = 0;

    virtual Complex operator*(const Current &c) const = 0;

    virtual void Print() const = 0;

    inline void SetP(const ATOOLS::Vec4D &p) { m_p=p; }

    inline void SetId(const Int_Vector &id) { m_id=id;   }
    inline void SetKey(const size_t &key)   { m_key=key; }

    inline ct::type Type() const { return m_type; }

    inline ATOOLS::Flavour Flav() const { return m_fl; }

    inline ATOOLS::Vec4D P() const { return m_p; }

    inline const Int_Vector &Id() const   { return m_id;   }
    inline size_t            Key() const  { return m_key;  }

    inline void AttachIn(Vertex *const v)  { m_in.push_back(v);  }
    inline void AttachOut(Vertex *const v) { m_out.push_back(v); }

  };// end of class Current

  std::ostream &operator<<(std::ostream &str,const Current &c);

}// end of namespace EXTRAXS

#endif
