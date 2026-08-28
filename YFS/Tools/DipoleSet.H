#ifndef YFS_Tools_DipoleSet_H
#define YFS_Tools_DipoleSet_H

#include "YFS/Tools/Dipole.H"

#include <functional>
#include <vector>
#include <cstddef>

namespace YFS {

  // Range over a subset of a DipoleSet, addressed by index so it stays valid
  // as the underlying vector grows.
  class DipoleView {
  public:
    class iterator {
    public:
      iterator(std::vector<Dipole> *d, const std::vector<std::size_t> *i,
               std::size_t n) : p_d(d), p_i(i), m_n(n) {}
      Dipole &operator*()  const { return (*p_d)[(*p_i)[m_n]]; }
      Dipole *operator->() const { return &(*p_d)[(*p_i)[m_n]]; }
      iterator &operator++() { ++m_n; return *this; }
      bool operator!=(const iterator &o) const { return m_n != o.m_n; }
      bool operator==(const iterator &o) const { return m_n == o.m_n; }
    private:
      std::vector<Dipole>            *p_d;
      const std::vector<std::size_t> *p_i;
      std::size_t                     m_n;
    };

    DipoleView(std::vector<Dipole> &d, const std::vector<std::size_t> &i)
      : p_d(&d), p_i(&i) {}

    iterator begin() const { return iterator(p_d, p_i, 0); }
    iterator end()   const { return iterator(p_d, p_i, p_i->size()); }
    std::size_t size()  const { return p_i->size(); }
    bool        empty() const { return p_i->empty(); }
    Dipole &operator[](std::size_t n) const { return (*p_d)[(*p_i)[n]]; }

  private:
    std::vector<Dipole>            *p_d;
    const std::vector<std::size_t> *p_i;
  };

  class DipoleSet {
  public:

    typedef std::function<double(const Dipole &)> ResonanceScore;

    // A dipole's born momenta are the momenta BEFORE that dipole radiates.
    // For the initial-initial dipole that is the born momenta of the process;
    // for the final-state dipoles it is the momenta after ISR but before FSR.
    // The two are therefore built at different points in the event, which is
    // why there are two entry points and not one.

    // Before ISR. Builds the initial-initial dipole.
    void BuildInitial(const ATOOLS::Flavour_Vector &flavs,
                      const ATOOLS::Vec4D_Vector   &momenta,
                      const ATOOLS::Vec4D_Vector   &born,
                      double alpha);

    // After ISR, before FSR. Builds the final-final and initial-final dipoles
    // with born = momenta, and selects the radiating FF subset. Leaves any
    // initial-initial dipole from BuildInitial untouched.
    void BuildFinal(const ATOOLS::Flavour_Vector &flavs,
                    const ATOOLS::Vec4D_Vector   &momenta,
                    double alpha, bool withIF,
                    const ResonanceScore &score = ResonanceScore());

    void Clear();

    // -- access -------------------------------------------------------------
    std::size_t size() const { return m_dipoles.size(); }
    bool empty()       const { return m_dipoles.empty(); }

    Dipole       &operator[](std::size_t i)       { return m_dipoles[i]; }
    const Dipole &operator[](std::size_t i) const { return m_dipoles[i]; }

    std::vector<Dipole>       &All()       { return m_dipoles; }
    const std::vector<Dipole> &All() const { return m_dipoles; }

    const std::vector<std::size_t> &Index(dipoletype::code t) const;

    const std::vector<std::size_t> &RadiatingIndex() const { return m_idxRad; }

    DipoleView ByType(dipoletype::code t) { return DipoleView(m_dipoles, Index(t)); }
    DipoleView FF() { return DipoleView(m_dipoles, m_idxFF); }
    DipoleView IF() { return DipoleView(m_dipoles, m_idxIF); }
    DipoleView Radiating() { return DipoleView(m_dipoles, m_idxRad); }

    Dipole       *II();
    const Dipole *II() const;

  private:
    // One charged leg, with its position in the event's flavour/momentum vectors.
    struct Leg {
      ATOOLS::Flavour flav;
      ATOOLS::Vec4D   mom, born;
      std::size_t     idx;
    };

    static std::vector<Leg> ChargedLegs(const ATOOLS::Flavour_Vector &flavs,
                                        const ATOOLS::Vec4D_Vector &momenta,
                                        const ATOOLS::Vec4D_Vector &born,
                                        std::size_t first, std::size_t last);
    void Add(const Leg &a, const Leg &b, dipoletype::code t, double alpha);
    void SelectRadiating(const ResonanceScore &score);

    std::vector<Dipole>      m_dipoles;
    std::vector<std::size_t> m_idxII, m_idxFF, m_idxIF, m_idxRad;
  };

}// end of namespace YFS

#endif
