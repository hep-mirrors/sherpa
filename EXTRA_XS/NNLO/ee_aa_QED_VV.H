#ifndef EXTRAXS__NNLO__ee_aa_QED_VV_H
#define EXTRAXS__NNLO__ee_aa_QED_VV_H
/*!
  \file ee_aa_QED_VV.H

  VV_Generator: Internal - the two-loop QED virtual of e+e- -> gamma gamma
  (hep-ph/0201274, EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H) as a loop provider
  for the double virtual (nlo_type vv, Process_Info::m_vvgenerator).

  Convention. Massless electrons, CDR, Catani's I^(1), I^(2) subtracted at
  mu^2 = the renormalisation scale set with SetRenScale (default s), so the
  result is finite by construction and no poles are reported
  (ProvidesPoles() false, IR() = IR2() = 0). In the provider normalisation
  of EXTRA_XS (DY_QCD_Virtual): the virtual is (alpha/2pi) x m_res x Born,
  here the two-loop finite remainder per Born in units of (alpha/2pi)^2,
    m_res.Finite() = [Finite^{0x2} + Finite^{1x1}]/T0 ,
  photonic only (no fermion loops, no light-by-light); multiply by alpha/2pi
  once more for the absolute size. ME_Born() is the spin-averaged tree
  2 e^4 (t/u + u/t) at the model's alpha_QED.

  A massive, YFS-subtracted use (ACRAIC/YFS BETA0_2) needs the pieces at
  mu^2 = s and its own massification; they are exposed, for the last Calc
  point, by the accessors below (YFS::AddPhotonicTwoLoopDiphoton reads them).
*/
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include <complex>

namespace EXTRAXS {

  class ee_aa_QED_VV : public PHASIC::Virtual_ME2_Base {
  private:
    //! last Calc point: x = -t/s, and the mu^2 = s pieces
    double m_x = 0.5, m_s = 0., m_fin02 = 0., m_fin11 = 0.;
    std::complex<double> m_r1;
    bool m_ok = false;
  public:
    ee_aa_QED_VV(const PHASIC::Process_Info &pi,
                 const ATOOLS::Flavour_Vector &flavs);

    void Calc(const ATOOLS::Vec4D_Vector &p) override;
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector &mom) override
    { return 4.*M_PI; }

    //! the last Calc point is e+e- -> gamma gamma with s > 0
    inline bool   Valid() const { return m_ok; }
    //! x = -t/s of the last point (symmetric under t <-> u)
    inline double X() const     { return m_x; }
    //! Finite^{0x2}/T0 (photonic two-loop interference), mu^2 = s
    inline double Fin02() const { return m_fin02; }
    //! Finite^{1x1}/T0 (one loop squared, photonic), mu^2 = s
    inline double Fin11() const { return m_fin11; }
    //! Catani-finite one-loop remainder, spin averaged, alpha/2pi, mu^2 = s
    inline std::complex<double> R1() const { return m_r1; }
    //! s of the last point
    inline double S() const     { return m_s; }
    /*!
      The last point's Laurent series (ee_aa_QED_Poles.H), per D-dimensional
      Born, (alpha/2pi)^L, D = 4 - 2 eps, mu^2 = s.
      Massless (MSbar normalisation): n02 = (P02 + Fin02)/T and
      n11 = (P11 + Fin11)/T, eps^-4 .. eps^0; m1 = 2 Re <M0|M1>/T,
      eps^-2 .. eps^2.
      Massive electrons of mass me (Mitov-Moch massification; OpenLoops/Recola
      normalisation (4pi)^eps/Gamma(1-eps), as Sherpa's YFS B): m1m eps^-1 ..
      eps^1, m2m eps^-2 .. eps^0. Only soft poles remain; YFS's subtraction
      cancels them as m1m_-1 = sigma_-1, m2m_-2 = sigma_-1^2/2,
      m2m_-1 = sigma_-1 m1m_0 (sigma = 2 B, alpha/2pi units).
    */
    void MasslessSeries(double n02[5], double n11[5], double m1[5]) const;
    void MassiveSeries(double me, double m1m[3], double m2m[3]) const;
  };

}

#endif
