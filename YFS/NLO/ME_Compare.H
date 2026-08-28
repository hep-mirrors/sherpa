#ifndef YFS_NLO_ME_Compare_H
#define YFS_NLO_ME_Compare_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

#include <map>
#include <string>
#include <fstream>

namespace ATOOLS {
class Histogram;
}

namespace YFS {

typedef std::map<std::string, ATOOLS::Histogram *> histograms1d;

// Point-by-point comparison of two matrix-element evaluations (either
// internal-vs-external or external-vs-external) for a real-emission process.
// Tracks mismatch statistics and histograms photon kinematics (energy, pT,
// angle to the nearest charged leg, angle to the beam axis, and - when
// there are 2+ photons, e.g. for RealReal - the minimum photon-photon
// opening angle) split by whether each point agreed with tolerance or not,
// to help correlate mismatches with the soft/collinear region, including
// the photon-photon collinear/soft region that only RealReal can reach.
// Shared between Real and RealReal so the diagnostic code isn't duplicated
// per class.
class ME_Compare {
public:
  ME_Compare(bool check, const std::string &histdir,
             const std::string &label, double tol = 1e-4);
  ~ME_Compare();

  // test: the value under test (e.g. internal ME, or a second external ME)
  // reference: the value it's being checked against (e.g. an external ME)
  // flavs: flavours of ALL legs (initial+final), same order as p
  // nin: number of initial-state legs in flavs/p
  void CheckAgreement(const ATOOLS::Vec4D_Vector &p, double test,
                      double reference, const ATOOLS::Flavour_Vector &flavs,
                      size_t nin);

private:
  bool m_check;
  double m_tol;
  std::string m_histdir, m_label;
  long int m_npoints, m_nbad, m_nbad_is, m_nbad_fs;
  double m_maxdev;
  histograms1d m_histograms1d;
  // per-photon scatter data (energy vs angle to nearest charged particle) for
  // spotting collinear enhancement of energetic photons in the large
  // corrections. All "bad" points are written; "ok" points capped to
  // m_scatter_max rows to bound the file size.
  std::ofstream m_scatter;
  long int m_scatter_ok, m_scatter_max;
};

}

#endif
