#ifndef YFS_Tools_DipoleSet_H
#define YFS_Tools_DipoleSet_H

#include "YFS/Main/Dipole.H"

#include <functional>
#include <memory>
#include <vector>
#include <cstddef>

namespace YFS {

  // Range over a subset of a DipoleSet, addressed by index so it stays valid
  // as the underlying vector grows.
  class DipoleView {
  public:
    class iterator {
    public:
      iterator(std::vector<std::unique_ptr<Dipole> > *d,
               const std::vector<std::size_t> *i,
               std::size_t n) : p_d(d), p_i(i), m_n(n) {}
      Dipole &operator*()  const { return *(*p_d)[(*p_i)[m_n]]; }
      Dipole *operator->() const { return (*p_d)[(*p_i)[m_n]].get(); }
      iterator &operator++() { ++m_n; return *this; }
      bool operator!=(const iterator &o) const { return m_n != o.m_n; }
      bool operator==(const iterator &o) const { return m_n == o.m_n; }
    private:
      std::vector<std::unique_ptr<Dipole> > *p_d;
      const std::vector<std::size_t>       *p_i;
      std::size_t                           m_n;
    };

    DipoleView(std::vector<std::unique_ptr<Dipole> > &d,
               const std::vector<std::size_t> &i) : p_d(&d), p_i(&i) {}

    iterator begin() const { return iterator(p_d, p_i, 0); }
    iterator end()   const { return iterator(p_d, p_i, p_i->size()); }
    std::size_t size()  const { return p_i->size(); }
    bool        empty() const { return p_i->empty(); }
    Dipole &operator[](std::size_t n) const { return *(*p_d)[(*p_i)[n]]; }

  private:
    std::vector<std::unique_ptr<Dipole> > *p_d;
    const std::vector<std::size_t>       *p_i;
  };

  class DipoleSet {
  public:

    typedef std::function<double(const Dipole &)> ResonanceScore;

    // A dipole's born momenta are the momenta BEFORE it radiates: the process
    // born momenta for II, the post-ISR pre-FSR momenta for the final-state
    // dipoles. Hence two entry points.

    // Before ISR. Builds the initial-initial dipole.
    void BuildInitial(const ATOOLS::Flavour_Vector &flavs,
                      const ATOOLS::Vec4D_Vector   &momenta,
                      const ATOOLS::Vec4D_Vector   &born,
                      double alpha);

    // After ISR, before FSR. Builds the final-final dipoles and selects the
    // radiating subset; leaves the II dipole untouched. born is the momenta
    // before FSR (usually the same vector as momenta).
    void BuildFinal(const ATOOLS::Flavour_Vector &flavs,
                    const ATOOLS::Vec4D_Vector   &momenta,
                    const ATOOLS::Vec4D_Vector   &born,
                    double alpha,
                    const ResonanceScore &score = ResonanceScore());

    // -- pole expansion (WW) ------------------------------------------------
    // The 4f final state of e+e- -> W+W- -> 4f resolved into the two W's and
    // their charged leptons, pairing each lepton with its partner neutrino
    // (|kf| differing by one, opposite sign).
    struct WWLegs {
      bool ok = false;                  //!< false => not a W+W- -> 4f state
      ATOOLS::Vec4D wm, wp;             //!< W-, W+ four-momenta (off-shell)
      std::size_t lm = 0, lp = 0;       //!< positions of the two charged leptons
      std::size_t nm = 0, np = 0;       //!< positions of the two neutrinos
    };

    static WWLegs FindWW(const ATOOLS::Flavour_Vector &flavs,
                         const ATOOLS::Vec4D_Vector   &momenta);

    bool BuildPole(const ATOOLS::Flavour_Vector &flavs,
                   const ATOOLS::Vec4D_Vector   &momenta,
                   const ATOOLS::Vec4D_Vector   &born,
                   double alpha, bool onshell, double maxdist,
                   bool emission);

    DipoleView Decay() { return DipoleView(m_dipoles, m_idxDec); }

    const WWLegs &WW() const { return m_ww; }

    // Initial-final dipoles: initial legs are the BORN beams, final legs are
    // at m_plab. The caller supplies the combined vector.
    void BuildIF(const ATOOLS::Flavour_Vector &flavs,
                 const ATOOLS::Vec4D_Vector   &momenta,
                 const ATOOLS::Vec4D_Vector   &born,
                 double alpha);

    void Clear();

    // -- access -------------------------------------------------------------
    std::size_t size() const { return m_dipoles.size(); }
    bool empty()       const { return m_dipoles.empty(); }

    Dipole       &operator[](std::size_t i)       { return *m_dipoles[i]; }
    const Dipole &operator[](std::size_t i) const { return *m_dipoles[i]; }


    const std::vector<std::size_t> &Index(dipoletype::code t) const;

    DipoleView ByType(dipoletype::code t) { return DipoleView(m_dipoles, Index(t)); }
    DipoleView FF() { return DipoleView(m_dipoles, m_idxFF); }
    DipoleView IF() { return DipoleView(m_dipoles, m_idxIF); }

    // At most one initial-initial dipole; II() throws if absent.
    bool HasII() const { return !m_idxII.empty(); }
    Dipole       &II();
    const Dipole &II() const;

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
    // Drop the named subsets, keeping the others.
    void DropTypes(bool dropII, bool dropFF, bool dropIF);
    void Add(const Leg &a, const Leg &b, dipoletype::code t, double alpha);
    // Add() for a pair that may contain a W leg; see the definition.
    void AddW(const Leg &a, const Leg &b, dipoletype::code t, double alpha,
              bool onshell, double MW);
    static void ProjectOnShell(ATOOLS::Vec4D &wm, ATOOLS::Vec4D &wp, double MW);
    void SelectRadiating(const ResonanceScore &score);

    // unique_ptr: a Dipole keeps its address (ISR caches a Dipole*).
    std::vector<std::unique_ptr<Dipole> > m_dipoles;
    std::vector<std::size_t> m_idxII, m_idxFF, m_idxIF, m_idxRad, m_idxDec;
    WWLegs m_ww;
  };

}// end of namespace YFS

#endif
