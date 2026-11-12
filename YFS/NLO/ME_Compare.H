#ifndef YFS_NLO_ME_Compare_H
#define YFS_NLO_ME_Compare_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

#include <map>
#include <string>
#include <fstream>
#include <memory>

namespace ATOOLS {
class Histogram;
}

namespace YFS {

typedef std::map<std::string, std::unique_ptr<ATOOLS::Histogram> > histograms1d;

//! Point-by-point comparison of two matrix-element evaluations of a
//! real-emission process. Tracks mismatch statistics and histograms photon
//! kinematics (energy, pT, angles, photon-photon opening angle) split by
//! agreement within tolerance.
class ME_Compare {
public:
  ME_Compare(bool check, const std::string &histdir,
             const std::string &label, double tol = 1e-4);
  ~ME_Compare();

  //! test: value under test, reference: value it is checked against,
  //! flavs: flavours of all legs in the order of p, nin: number of initial legs
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
  //! per-photon scatter data (energy vs angle to nearest charged particle);
  //! all bad points are written, ok points capped at m_scatter_max rows
  std::ofstream m_scatter;
  long int m_scatter_ok, m_scatter_max;
};

}

#endif
