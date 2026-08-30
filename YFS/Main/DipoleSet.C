#include "YFS/Main/DipoleSet.H"

#include "ATOOLS/Math/Poincare.H"
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
    m_idxDec.clear();
  }

  const std::vector<std::size_t> &DipoleSet::Index(dipoletype::code t) const {
    switch (t) {
    case dipoletype::initial: return m_idxII;
    case dipoletype::final:   return m_idxFF;
    case dipoletype::ifi:     return m_idxIF;
    }
    return m_idxII;
  }

  Dipole &DipoleSet::II() {
    if (m_idxII.empty())
      THROW(fatal_error, "No initial-initial dipole; DipoleSet::BuildInitial "
                         "has not run, or the initial state carries no charge.");
    return *m_dipoles[m_idxII.front()];
  }

  const Dipole &DipoleSet::II() const {
    if (m_idxII.empty())
      THROW(fatal_error, "No initial-initial dipole; DipoleSet::BuildInitial "
                         "has not run, or the initial state carries no charge.");
    return *m_dipoles[m_idxII.front()];
  }

  void DipoleSet::DropTypes(bool dropII, bool dropFF, bool dropIF) {
    // Move the survivors, so every kept Dipole keeps its address.
    std::vector<std::unique_ptr<Dipole> > keep;
    std::vector<dipoletype::code> kind;
    auto take=[&](const std::vector<std::size_t> &idx, bool drop, dipoletype::code t){
      if (drop) return;
      for (std::size_t i : idx) { keep.push_back(std::move(m_dipoles[i])); kind.push_back(t); }
    };
    // The decay dipoles have their own index list, so they are tracked
    // separately here and re-registered below. They are part of the pole
    // expansion's final-state structure, so they die with the final-final
    // dipoles: BuildFinal runs at the top of every event, and a decay dipole
    // that outlived its event would contribute a form-factor term built from
    // the previous event's momenta to an event that fell back to the flat
    // scheme.
    std::vector<std::unique_ptr<Dipole> > keepdec;
    if (!dropFF)
      for (std::size_t i : m_idxDec) keepdec.push_back(std::move(m_dipoles[i]));
    take(m_idxII, dropII, dipoletype::initial);
    take(m_idxFF, dropFF, dipoletype::final);
    take(m_idxIF, dropIF, dipoletype::ifi);
    Clear();
    for (std::size_t n(0); n < keep.size(); ++n) {
      m_dipoles.push_back(std::move(keep[n]));
      const std::size_t i(m_dipoles.size() - 1);
      switch (kind[n]) {
      case dipoletype::initial: m_idxII.push_back(i); break;
      case dipoletype::final:   m_idxFF.push_back(i);
                                if (m_dipoles[i]->IsResonance()) m_idxRad.push_back(i); break;
      case dipoletype::ifi:     m_idxIF.push_back(i); break;
      }
    }
    for (std::size_t n(0); n < keepdec.size(); ++n) {
      m_dipoles.push_back(std::move(keepdec[n]));
      m_idxDec.push_back(m_dipoles.size() - 1);
    }
  }

  void DipoleSet::Add(const Leg &a, const Leg &b, dipoletype::code t,
                      double alpha) {
    Flavour_Vector fl{a.flav, b.flav};
    Vec4D_Vector   mo{a.mom,  b.mom};
    Vec4D_Vector   bo{a.born, b.born};
    std::unique_ptr<Dipole> D(new Dipole(fl, mo, bo, t, alpha));
    D->SetFlavLab((int)a.idx, (int)b.idx);
    m_dipoles.push_back(std::move(D));
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

    DropTypes(true, false, false);

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
                             const Vec4D_Vector   &born,
                             double alpha,
                             const ResonanceScore &score) {
    if (momenta.size() != flavs.size() || born.size() != flavs.size())
      THROW(fatal_error, "Inconsistent vector sizes in DipoleSet::BuildFinal");

    DropTypes(false, true, false);

    const std::vector<Leg> in (ChargedLegs(flavs, momenta, born, 0, 2));
    const std::vector<Leg> out(ChargedLegs(flavs, momenta, born, 2, flavs.size()));

    // Every unique charged pair: the virtual and form-factor sums need all of
    // them. Which ones radiate is decided by SelectRadiating.
    for (std::size_t i(0); i < out.size(); ++i)
      for (std::size_t j(i + 1); j < out.size(); ++j)
        Add(out[i], out[j], dipoletype::final, alpha);

    SelectRadiating(score);
  }

  void DipoleSet::BuildIF(const Flavour_Vector &flavs,
                          const Vec4D_Vector   &momenta,
                          const Vec4D_Vector   &born,
                          double alpha) {
    if (momenta.size() != flavs.size() || born.size() != flavs.size())
      THROW(fatal_error, "Inconsistent vector sizes in DipoleSet::BuildIF");

    DropTypes(false, false, true);

    const std::vector<Leg> in (ChargedLegs(flavs, momenta, born, 0, 2));
    const std::vector<Leg> out(ChargedLegs(flavs, momenta, born, 2, flavs.size()));

    for (std::size_t i(0); i < in.size(); ++i)
      for (std::size_t j(0); j < out.size(); ++j)
        Add(in[i], out[j], dipoletype::ifi, alpha);
  }

  void DipoleSet::SelectRadiating(const ResonanceScore &score) {
    m_idxRad.clear();
    for (std::size_t i : m_idxFF) m_dipoles[i]->SetResonance(false);

    std::set<int> used;
    for (int pass(0); pass < 2; ++pass) {
      std::vector<std::pair<double, std::size_t> > cand;
      for (std::size_t i : m_idxFF) {
        Dipole &D(*m_dipoles[i]);
        if (D.m_QiQj >= 0) continue;                    // opposite charges only
        if (pass == 0 && !D.IsDecayAllowed()) continue;  // same flavour first
        if (used.count(D.Left()) || used.count(D.Right())) continue;
        const double s(score ? score(D) : std::numeric_limits<double>::max());
        cand.push_back(std::make_pair(s, i));
      }
      std::stable_sort(cand.begin(), cand.end());
      for (std::size_t c(0); c < cand.size(); ++c) {
        Dipole &D(*m_dipoles[cand[c].second]);
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
      legs.insert(m_dipoles[i]->Left());
      legs.insert(m_dipoles[i]->Right());
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

  DipoleSet::WWLegs DipoleSet::FindWW(const Flavour_Vector &flavs,
                                      const Vec4D_Vector   &momenta) {
    WWLegs w;
    // Pair each charged lepton with the neutrino of its own generation and
    // opposite sign: l-(kf=+c) goes with nubar_l(kf=-(c+1)), l+(kf=-c) with
    // nu_l(kf=+(c+1)). That is the W it came from, whatever the other W did.
    for (std::size_t i(2); i < flavs.size(); ++i) {
      const long kf(flavs[i].Kfcode());
      if (kf != 11 && kf != 13 && kf != 15) continue;
      const bool neg(!flavs[i].IsAnti());          // l- rather than l+
      const long nukf(kf + 1);
      for (std::size_t j(2); j < flavs.size(); ++j) {
        if (j == i || flavs[j].Kfcode() != nukf) continue;
        if (flavs[j].IsAnti() != neg) continue;    // l- with nubar, l+ with nu
        if (neg) { w.wm = momenta[i] + momenta[j]; w.lm = i; w.nm = j; }
        else     { w.wp = momenta[i] + momenta[j]; w.lp = i; w.np = j; }
        break;
      }
    }
    // Both W's must have been found, and between them they must account for
    // the whole final state -- otherwise this is not W+W- -> 4f and the pole
    // expansion does not apply.
    w.ok = (w.lm != w.lp) && !IsZero(w.wm.E()) && !IsZero(w.wp.E()) &&
           flavs.size() == 6;
    return w;
  }

  bool DipoleSet::BuildPole(const Flavour_Vector &flavs,
                            const Vec4D_Vector   &momenta,
                            const Vec4D_Vector   &born,
                            double alpha, bool onshell, double maxdist,
                            bool emission) {
    if (momenta.size() != flavs.size() || born.size() != flavs.size())
      THROW(fatal_error, "Inconsistent vector sizes in DipoleSet::BuildPole");

    const WWLegs w(FindWW(flavs, momenta));
    if (!w.ok) return false;
    const WWLegs wb(FindWW(flavs, born));
    if (!wb.ok) return false;

    // Leading-pole validity. With no generation cuts the integrator visits
    // configurations where the "W" reconstructs to a few GeV; calling that a
    // narrow resonance and exponentiating radiation off it is meaningless, and
    // it is where the decay-stage log(M/m_l) goes wild. Both W's must be near
    // the pole or the whole event goes back to the flat scheme.
    const double MW0(Flavour(kf_Wplus).Mass()), GW(Flavour(kf_Wplus).Width());
    if (GW > 0.) {
      const double dm(fabs(w.wm.Mass() - MW0)/GW), dp(fabs(w.wp.Mass() - MW0)/GW);
      if (dm > maxdist || dp > maxdist) return false;
    }

    // Drop whatever flat final-state structure is there; keep the II dipole,
    // which is the same (e-,e+) pair in both schemes.
    DropTypes(false, true, true);
    m_idxDec.clear();

    const Flavour Wm(Flavour(kf_Wplus).Bar()), Wp(kf_Wplus);
    const double MW(Flavour(kf_Wplus).Mass());

    // On-shell variation: rescale the spatial momentum in the W-pair rest
    // frame so each W sits exactly at M_W. The nominal leaves them off shell,
    // carrying the invariant mass of their own decay products.
    Vec4D wmm(w.wm), wpm(w.wp), wmb(wb.wm), wpb(wb.wp);
    if (onshell) {
      ProjectOnShell(wmm, wpm, MW);
      ProjectOnShell(wmb, wpb, MW);
    }

    // Leg indices are reported as the positions of the charged leptons the W's
    // decayed to: the W is not itself an entry in the event's flavour vector,
    // and this is what makes the dump readable.
    m_ww = w;

    const Leg lwm{Wm, wmm, wmb, w.lm};
    const Leg lwp{Wp, wpm, wpb, w.lp};

    const std::vector<Leg> in(ChargedLegs(flavs, momenta, born, 0, 2));

    // -- production stage --------------------------------------------------
    AddW(lwm, lwp, dipoletype::final, alpha, onshell, MW);
    for (std::size_t i(0); i < in.size(); ++i) {
      AddW(in[i], lwm, dipoletype::ifi, alpha, onshell, MW);
      AddW(in[i], lwp, dipoletype::ifi, alpha, onshell, MW);
    }
    // The W pair is the only thing that radiates at the production stage.
    // With emission off nothing in the final state radiates and the pole
    // scheme changes the exponent alone.
    for (std::size_t i : m_idxFF) m_dipoles[i]->SetResonance(emission);
    m_idxRad.clear();
    if (emission) m_idxRad = m_idxFF;

    // -- decay stage -------------------------------------------------------
    // One dipole per W: (W, charged lepton). The neutrino is neutral, so each
    // leptonic decay is a two-leg dipole and nothing else. Built as ifi so the
    // parent W carries the opposite theta to its daughter, which is what makes
    // the pair charge conserving (Q_W = Q_l for a leptonic decay, with the W
    // incoming) and reproduces eq (12) of hep-ph/0302065.
    // The lepton legs keep their own event positions; the W legs above reuse
    // them, so a decay dipole prints as legs=(l,l). Offset the W's label by the
    // event size so the dump distinguishes parent from daughter.
    const Leg dlm{flavs[w.lm], momenta[w.lm], born[w.lm], w.lm};
    const Leg dlp{flavs[w.lp], momenta[w.lp], born[w.lp], w.lp};
    const Leg pwm{Wm, wmm, wmb, w.lm + flavs.size()};
    const Leg pwp{Wp, wpm, wpb, w.lp + flavs.size()};
    const std::size_t first(m_dipoles.size());
    AddW(pwm, dlm, dipoletype::ifi, alpha, onshell, MW);
    AddW(pwp, dlp, dipoletype::ifi, alpha, onshell, MW);
    // Move the two just-added dipoles out of the interference bucket: they are
    // a separate stage, not production-stage IF pairs, and summing them
    // together would double count.
    for (std::size_t i(first); i < m_dipoles.size(); ++i) {
      m_idxIF.erase(std::remove(m_idxIF.begin(), m_idxIF.end(), i),
                    m_idxIF.end());
      m_idxDec.push_back(i);
    }
    return true;
  }

  // Add(), then fix up any W leg's mass to the invariant mass of the momentum
  // it was actually built from. The Dipole constructor takes masses from the
  // Flavour, i.e. on shell; off shell that would leave GetMass() and
  // p.Mass() disagreeing for the same leg, and the form factors use both.
  void DipoleSet::AddW(const Leg &a, const Leg &b, dipoletype::code t,
                       double alpha, bool onshell, double MW) {
    Add(a, b, t, alpha);
    if (onshell) return;
    Dipole &D(*m_dipoles.back());
    for (int i(0); i < 2; ++i) {
      if (!D.GetFlav(i).IsVector()) continue;
      const double m2(D.GetBornMomenta(i).Abs2());
      if (m2 > 0.) D.SetMass(i, sqrt(m2));
    }
  }

  // Rescale both W momenta in their pair rest frame so each sits at M_W,
  // keeping the pair total and the directions. Used only for the on-shell
  // variation.
  void DipoleSet::ProjectOnShell(Vec4D &wm, Vec4D &wp, double MW) {
    const Vec4D Q(wm + wp);
    const double s(Q.Abs2());
    if (s <= 4.*MW*MW) return;            // below the on-shell pair threshold
    Poincare boost(Q);
    Vec4D a(wm), b(wp);
    boost.Boost(a); boost.Boost(b);
    const double p(0.5*sqrt(s - 4.*MW*MW));
    const Vec3D n(Vec3D(a)/Vec3D(a).Abs());
    a = Vec4D(sqrt(s)/2., p*n);
    b = Vec4D(sqrt(s)/2., -p*n);
    boost.BoostBack(a); boost.BoostBack(b);
    wm = a; wp = b;
  }

}// end of namespace YFS
