#ifndef ATOOLS_Phys_SelectorList
#define ATOOLS_Phys_SelectorList

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"

#include <vector>

namespace ATOOLS {

  class Selector_Particle {
  private:
    ATOOLS::Flavour m_fl;
    ATOOLS::Vec4D   m_p;
  public:
    Selector_Particle() :
      m_fl(ATOOLS::Flavour(kf_none)), m_p(ATOOLS::Vec4D(0.,0.,0.,0.)) {}
    Selector_Particle(ATOOLS::Flavour fl, ATOOLS::Vec4D p) :
      m_fl(fl), m_p(p) {}
    ~Selector_Particle();

    inline ATOOLS::Flavour Flavour()  const { return m_fl; }
    inline ATOOLS::Vec4D   Momentum() const { return m_p; }

    inline void SetFlavour(ATOOLS::Flavour fl) { m_fl=fl; }
    inline void SetMomentum(ATOOLS::Vec4D p)   { m_p=p; }
  };

  std::ostream &operator<<(std::ostream &ostr,const Selector_Particle &p);

  class Selector_List : public std::vector<Selector_Particle> {
  private:
    size_t m_nin;
  public:
    Selector_List();
    Selector_List(const ATOOLS::Flavour_Vector &fl,
                  const ATOOLS::Vec4D_Vector &p,
                  const size_t& nin);
    Selector_List(const ATOOLS::Flavour *fl, size_t n,
                  const ATOOLS::Vec4D_Vector &p,
                  const size_t& nin);

    void AddFlavour(ATOOLS::Flavour fl);

    void SetMomenta(ATOOLS::Vec4D_Vector p);
    ATOOLS::Vec4D_Vector ExtractMomenta();

    inline void SetNIn(const size_t& nin) { m_nin=nin; }
    inline const size_t& NIn() const      { return m_nin; }
  };

  std::ostream &operator<<(std::ostream &ostr,const Selector_List &p);
}

#endif
