#include "YFS/Main/YFS_FB_Analysis.H"

#include "ATOOLS/Math/Poincare.H"

#include <cmath>
#include <utility>

using namespace ATOOLS;

namespace YFS {

  std::string FBTag(fbdef::code def) {
    switch (def) {
      case fbdef::lab:           return "Lab";
      case fbdef::pair_rest:     return "PR";
      case fbdef::collins_soper: return "CS";
    }
    return "Unknown";
  }

  YFS_FB_Analysis::YFS_FB_Analysis(const std::vector<fbdef::code>& enabled,
                                   long int target_kf)
    : m_defs(enabled), m_target_kf(std::labs(target_kf)) {
    if (m_defs.empty())
      m_defs = {fbdef::lab, fbdef::pair_rest, fbdef::collins_soper};
  }

  void YFS_FB_Analysis::ReferenceIndices(const Flavour_Vector& fl,
                                         int& neg, int& pos) const {
    neg = pos = -1;
    for (size_t i = 2; i < fl.size(); ++i) {
      if (fl[i].IsPhoton()) continue;
      if (m_target_kf && (long int)fl[i].Kfcode() != m_target_kf) continue;
      const double q = fl[i].Charge();
      if (q < 0.0 && neg < 0) neg = (int)i;
      if (q > 0.0 && pos < 0) pos = (int)i;
    }
  }

  Vec4D YFS_FB_Analysis::IncomingMinus(const Vec4D_Vector& p,
                                       const Flavour_Vector& fl) const {
    for (size_t i = 0; i < 2 && i < fl.size(); ++i)
      if (fl[i].Charge() < 0.0) return p[i];
    if (p.size() >= 2) return (p[0][3] >= 0.0) ? p[0] : p[1];
    return Vec4D(1.0, 0.0, 0.0, 1.0);
  }

  double YFS_FB_Analysis::CosTheta(fbdef::code def,
                                   const Vec4D_Vector& p,
                                   const Flavour_Vector& fl) const {
    int neg, pos;
    ReferenceIndices(fl, neg, pos);
    if (neg < 0) return -2.0;

    Vec4D ref  = p[neg];
    Vec4D emin = IncomingMinus(p, fl);

    if (def == fbdef::lab)
      return ref.CosTheta(emin);

    // Boost the reference particle (and beams) into the outgoing pair rest
    // frame. The pair is the reference particle plus its positive partner if
    // present, otherwise all non-photon final-state momenta.
    Vec4D Q;
    if (pos >= 0) {
      Q = p[neg] + p[pos];
    } else {
      for (size_t i = 2; i < p.size(); ++i)
        if (!fl[i].IsPhoton()) Q += p[i];
    }
    if (Q.PSpat2() < 0.0 || Q.Abs2() <= 0.0) return -2.0;

    Poincare boost(Q);
    boost.Boost(ref);

    if (def == fbdef::pair_rest) {
      Vec4D emin_rest = emin;
      boost.Boost(emin_rest);
      return ref.CosTheta(emin_rest);
    }

    // Collins-Soper: axis bisecting the two beam directions in the pair rest
    // frame. With ISR the beams are not back-to-back there, so the bisector
    // differs from either beam direction.
    if (def == fbdef::collins_soper) {
      if (p.size() < 2) return -2.0;
      Vec4D b1 = p[0], b2 = p[1];
      boost.Boost(b1);
      boost.Boost(b2);
      const double n1 = b1.PSpat();
      const double n2 = b2.PSpat();
      if (n1 <= 0.0 || n2 <= 0.0) return -2.0;
      // Orient beam 1 along the incoming e- so the sign convention matches.
      const double s1 = (fl.size() > 0 && fl[0].Charge() < 0.0) ? 1.0 : -1.0;
      const double s2 = -s1;
      Vec4D axis(0.0,
                 s1 * b1[1] / n1 + s2 * b2[1] / n2,
                 s1 * b1[2] / n1 + s2 * b2[2] / n2,
                 s1 * b1[3] / n1 + s2 * b2[3] / n2);
      if (axis.PSpat2() <= 0.0) return -2.0;
      return ref.CosTheta(axis);
    }

    return -2.0;
  }

  void YFS_FB_Analysis::SplitWeights(Weights& w, const Vec4D_Vector& p,
                                     const Flavour_Vector& fl) const {
    // Snapshot the existing named contributions before adding new ones, since
    // Weights::operator[] appends. Index 0 is the (unnamed) nominal entry.
    std::vector<std::pair<std::string, double>> base;
    for (size_t i = 1; i < w.Size(); ++i)
      base.emplace_back(w.Name(i), w[i]);
    if (base.empty()) return;

    for (auto def : m_defs) {
      const double cth = CosTheta(def, p, fl);
      const bool valid = (cth >= -1.0);
      const bool forward = valid && (cth > 0.0);
      const std::string tag = FBTag(def);
      for (const auto& c : base) {
        w[c.first + "_" + tag + "_F"] = (valid &&  forward) ? c.second : 0.0;
        w[c.first + "_" + tag + "_B"] = (valid && !forward) ? c.second : 0.0;
      }
    }
  }

}  // namespace YFS
