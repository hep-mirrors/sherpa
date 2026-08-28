#include "YFS/Tools/DipoleSet.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>
#include <limits>
#include <set>
#include <utility>

using namespace ATOOLS;

namespace YFS {

  std::vector<DipoleSet::Leg>
  DipoleSet::ChargedLegs(const Flavour_Vector &flavs,
                         const Vec4D_Vector &momenta,
                         const Vec4D_Vector &born,
                         std::size_t first, std::size_t last) {
    std::vector<Leg> legs;
    for (std::size_t i(first); i < last; ++i) {
      if (flavs[i].IntCharge() == 0) continue;
      legs.push_back(Leg{flavs[i], momenta[i], born[i], i});
    }
    return legs;
  }

  void DipoleSet::Clear() {
    m_dipoles.clear();
    m_idxII.clear(); m_idxFF.clear(); m_idxIF.clear(); m_idxRad.clear();
  }

  const std::vector<std::size_t> &DipoleSet::Index(dipoletype::code t) const {
    switch (t) {
    case dipoletype::initial: return m_idxII;
    case dipoletype::final:   return m_idxFF;
    case dipoletype::ifi:     return m_idxIF;
    }
    return m_idxII;
  }

  Dipole *DipoleSet::II() {
    return m_idxII.empty() ? nullptr : &m_dipoles[m_idxII.front()];
  }

  const Dipole *DipoleSet::II() const {
    return m_idxII.empty() ? nullptr : &m_dipoles[m_idxII.front()];
  }

  void DipoleSet::Add(const Leg &a, const Leg &b, dipoletype::code t,
                      double alpha) {
    Flavour_Vector fl{a.flav, b.flav};
    Vec4D_Vector   mo{a.mom,  b.mom};
    Vec4D_Vector   bo{a.born, b.born};
    Dipole D(fl, mo, bo, t, alpha);
    D.SetFlavLab((int)a.idx, (int)b.idx);
    m_dipoles.push_back(D);
    const std::size_t i(m_dipoles.size() - 1);
    switch (t) {
    case dipoletype::initial: m_idxII.push_back(i); break;
    case dipoletype::final:   m_idxFF.push_back(i); break;
    case dipoletype::ifi:     m_idxIF.push_back(i); break;
    }
  }

  void DipoleSet::BuildInitial(const Flavour_Vector &flavs,
                               const Vec4D_Vector   &momenta,
                               const Vec4D_Vector   &born,
                               double alpha) {
    if (momenta.size() != flavs.size() || born.size() != flavs.size())
      THROW(fatal_error, "Inconsistent vector sizes in DipoleSet::BuildInitial");
    if (flavs.size() < 2)
      THROW(fatal_error, "DipoleSet::BuildInitial needs at least two particles");

    Clear();

    const std::vector<Leg> in (ChargedLegs(flavs, momenta, born, 0, 2));
    if (in.size() != 2) {
      if (!in.empty())
        msg_Error() << METHOD << "(): expected two charged initial-state "
                    << "particles, found " << in.size()
                    << "; no initial-initial dipole built." << std::endl;
      return;
    }
    Add(in[0], in[1], dipoletype::initial, alpha);
  }

  void DipoleSet::BuildFinal(const Flavour_Vector &flavs,
                             const Vec4D_Vector   &momenta,
                             double alpha, bool withIF,
                             const ResonanceScore &score) {
    if (momenta.size() != flavs.size())
      THROW(fatal_error, "Inconsistent vector sizes in DipoleSet::BuildFinal");

    // Keep the initial-initial dipole built before ISR; drop and rebuild
    // everything downstream of it.
    std::vector<Dipole> keep;
    for (std::size_t i : m_idxII) keep.push_back(m_dipoles[i]);
    Clear();
    for (const Dipole &D : keep) {
      m_dipoles.push_back(D);
      m_idxII.push_back(m_dipoles.size() - 1);
    }

    // born = momenta: these dipoles have not radiated yet.
    const std::vector<Leg> in (ChargedLegs(flavs, momenta, momenta, 0, 2));
    const std::vector<Leg> out(ChargedLegs(flavs, momenta, momenta, 2, flavs.size()));

    // Every unique charged pair: the virtual and form-factor sums need all of
    // them. Which ones radiate is decided by SelectRadiating.
    for (std::size_t i(0); i < out.size(); ++i)
      for (std::size_t j(i + 1); j < out.size(); ++j)
        Add(out[i], out[j], dipoletype::final, alpha);

    if (withIF)
      for (std::size_t i(0); i < in.size(); ++i)
        for (std::size_t j(0); j < out.size(); ++j)
          Add(in[i], out[j], dipoletype::ifi, alpha);

    SelectRadiating(score);
  }

  void DipoleSet::SelectRadiating(const ResonanceScore &score) {
    m_idxRad.clear();
    for (std::size_t i : m_idxFF) m_dipoles[i].SetResonance(false);

    std::set<int> used;
    for (int pass(0); pass < 2; ++pass) {
      std::vector<std::pair<double, std::size_t> > cand;
      for (std::size_t i : m_idxFF) {
        Dipole &D(m_dipoles[i]);
        if (D.m_QiQj >= 0) continue;                    // opposite charges only
        if (pass == 0 && !D.IsDecayAllowed()) continue;  // same flavour first
        if (used.count(D.Left()) || used.count(D.Right())) continue;
        const double s(score ? score(D) : std::numeric_limits<double>::max());
        cand.push_back(std::make_pair(s, i));
      }
      std::stable_sort(cand.begin(), cand.end());
      for (std::size_t c(0); c < cand.size(); ++c) {
        Dipole &D(m_dipoles[cand[c].second]);
        if (used.count(D.Left()) || used.count(D.Right())) continue;
        D.SetResonance(true);
        m_idxRad.push_back(cand[c].second);
        used.insert(D.Left());
        used.insert(D.Right());
        msg_Debugging() << METHOD << "(): radiating dipole (" << D.Left() << ","
                        << D.Right() << ") " << D.m_flavs[0] << D.m_flavs[1]
                        << " pass=" << pass << " score=" << cand[c].first
                        << std::endl;
      }
    }

    std::set<int> legs;
    for (std::size_t i : m_idxFF) {
      legs.insert(m_dipoles[i].Left());
      legs.insert(m_dipoles[i].Right());
    }
    for (int l : legs) {
      if (used.count(l)) continue;
      static bool warned(false);
      if (!warned) {
        msg_Error() << METHOD << "(): charged final-state particle at position "
                    << l << " enters no radiating dipole, so the final-state "
                    << "eikonal current does not conserve charge. This is "
                    << "expected for an odd number of charged final-state "
                    << "particles. Further warnings suppressed." << std::endl;
        warned = true;
      }
      msg_Debugging() << METHOD << "(): unpaired charged leg at " << l << std::endl;
    }
  }

}// end of namespace YFS
