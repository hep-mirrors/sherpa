#ifndef YFS_Main_Soft_Region_H
#define YFS_Main_Soft_Region_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/MathTools.H"
#include <string>
#include <sstream>

namespace YFS {

  /*!
    The soft region Omega of a YFS form factor or crude density: the photons
    with energy below Emin in one frame,
      Omega = { k : k.n < Emin } ,  n the frame's unit time direction.
    B-tilde(Omega) is a Lorentz scalar for a given Omega
    and what the frame changes is the region. The energy-based formula
    (YFS_Form_Factor::BVR_full(p1p2, E1, E2, ...)) takes the leg energies IN
    the region's frame, so the frame used to be implicit in where the caller
    took E1, E2 from. A SoftRegion names the frame
    and computes the energies itself.

    Three kinds, each reproducing the arithmetic of the code it replaced:
      Lab       the frame the momenta are given in (the beams' rest frame after
                YFS_Handler's CMS boost, "P"): E = p[0].
      RestOf    the rest frame of a timelike vector V: legs boosted with
                Poincare(V) (no boost when V is at rest to 1e-14).
      PairRest  the rest frame of a two-body system of invariant mass^2 s,
                energies from two-body kinematics,
                E1 = (s + m1^2 - m2^2)/(2 sqrt s) (FSR's Q frame, where the
                energies are computed, not carried by a vector).
  */
  struct SoftRegion {
    enum class kind { lab, rest_of, pair_rest };
    kind        m_kind = kind::lab;
    double      m_emin = 0.;
    ATOOLS::Vec4D m_frame;      //!< rest_of: the vector whose rest frame it is
    double      m_s = 0.;       //!< pair_rest: the pair's invariant mass^2
    bool        m_threshold_fallback = false;  //!< pair_rest: E = sqrt(s)/2 below threshold
    std::string m_tag;          //!< for the ledger (DumpDipoles)

    static SoftRegion Lab(double emin, const std::string &tag = "P") {
      SoftRegion r; r.m_kind = kind::lab; r.m_emin = emin; r.m_tag = tag; return r;
    }
    static SoftRegion RestOf(const ATOOLS::Vec4D &V, double emin,
                             const std::string &tag = "rest") {
      SoftRegion r; r.m_kind = kind::rest_of; r.m_frame = V; r.m_emin = emin;
      r.m_tag = tag; return r;
    }
    static SoftRegion PairRest(double s, double emin, bool threshold_fallback,
                               const std::string &tag = "Q") {
      SoftRegion r; r.m_kind = kind::pair_rest; r.m_s = s; r.m_emin = emin;
      r.m_threshold_fallback = threshold_fallback; r.m_tag = tag; return r;
    }

    double Emin() const { return m_emin; }

    /*!
      The legs in the region's frame for the form factor: their invariant
      product, and their energies there. p1p2 is taken from the vectors the
      formula used before (boosted ones for rest_of, the given ones
      otherwise), so that every caller stays bit-identical.
    */
    void Legs(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2, double m1, double m2,
              double &p1p2, double &E1, double &E2) const {
      switch (m_kind) {
      case kind::lab:
        p1p2 = p1*p2; E1 = p1[0]; E2 = p2[0]; return;
      case kind::rest_of: {
        ATOOLS::Vec4D r1(p1), r2(p2);
        if (m_frame.Abs2() > 0. && !ATOOLS::IsZero(m_frame.PSpat()/m_frame[0], 1e-14)) {
          ATOOLS::Poincare rest(m_frame);
          rest.Boost(r1);
          rest.Boost(r2);
        }
        p1p2 = r1*r2; E1 = r1.E(); E2 = r2.E(); return;
      }
      case kind::pair_rest:
        p1p2 = p1*p2;
        if (m_threshold_fallback && !(m_s > ATOOLS::sqr(m1 + m2))) {
          E1 = E2 = 0.5*sqrt(m_s); return;
        }
        E1 = (m_s + m1*m1 - m2*m2)/(2.*sqrt(m_s));
        E2 = (m_s + m2*m2 - m1*m1)/(2.*sqrt(m_s));
        return;
      }
    }

    //! a photon's energy in the region's frame (the ledger check)
    double PhotonEnergy(const ATOOLS::Vec4D &k) const {
      switch (m_kind) {
      case kind::lab: return k[0];
      case kind::rest_of: {
        ATOOLS::Vec4D r(k);
        if (m_frame.Abs2() > 0. && !ATOOLS::IsZero(m_frame.PSpat()/m_frame[0], 1e-14)) {
          ATOOLS::Poincare rest(m_frame); rest.Boost(r);
        }
        return r[0];
      }
      case kind::pair_rest: return -1.;   // no vector: not checkable per photon
      }
      return -1.;
    }

    std::string Describe() const {
      std::ostringstream o;
      o<<m_tag<<(m_kind == kind::lab ? "-sphere" : m_kind == kind::rest_of ? "-rest-sphere"
                 : "-pair-sphere")<<" Emin="<<m_emin;
      return o.str();
    }
  };

}

#endif
