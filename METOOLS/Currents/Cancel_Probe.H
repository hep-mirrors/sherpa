#ifndef METOOLS_Currents_Cancel_Probe_H
#define METOOLS_Currents_Cancel_Probe_H
// DIAGNOSTIC (Phase 0). Records the worst cancellation seen while summing
// current contributions: the ratio max(|a|,|b|)/|a+b|. The measured error law
// dev ~ eps*s/(2p.k) is what an absolute error eps*s in a quantity of true size
// 2p.k looks like, so a ratio of order s/(2p.k) ~ 1e13 here identifies the site.
#include <cstddef>
namespace METOOLS {
  struct Cancel_Probe {
    static double s_worst;      // worst ratio this event
    static double s_worst_mag;  // magnitude of the largest term involved
    static size_t s_n;          // sums examined
    static void Reset() { s_worst=0.0; s_worst_mag=0.0; s_n=0; }
    static void Note(const double &amax,const double &res) {
      ++s_n;
      if (res<=0.0 || amax<=0.0) return;
      const double r(amax/res);
      if (r>s_worst) { s_worst=r; s_worst_mag=amax; }
    }
  };
}
#endif
