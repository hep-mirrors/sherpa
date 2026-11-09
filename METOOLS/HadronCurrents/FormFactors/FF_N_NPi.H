#ifndef METOOLS_Current_Library_FF_N_NPi_H
#define METOOLS_Current_Library_FF_N_NPi_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  class FF_N_N;

  //////////////////////////////////////////////////////////////////////////
  /// Strong form factors for VA_N_NPi, evaluated at a virtuality x that the
  /// current supplies (x = -k^2 of the pion, or Q^2 / -t for FPiPole):
  ///   FpiNN   : pi-N-N coupling g_piNN(x), WITHOUT the isospin factor c_pi,
  ///             which depends on the vertex and is applied in VA_N_NPi.
  ///   FPiPole : form factor of the pion-pole part of the weak contact term.
  /// FORM_FACTOR:
  ///   FpiNN   : 0 zero, 1 point (Goldberger-Treiman g_A m_N/f_pi),
  ///             10 PionDipole,
  ///             20-25 PCAC, (m_N/f_pi) F_A(x) with F_A from FF_N_N (same
  ///             numbers as AXIAL_FF there).
  ///   FPiPole : 0 zero, 1 point (= 1),
  ///             11 rho dominance, m_rho^2/(m_rho^2 + x) (default),
  ///             12 F1V(x) = F1^p - F1^n from FF_N_N (VECTOR_FF model).
  /// The m_pi pair must be a proton and a neutron index, so that the
  /// internal FF_N_N ("FA" or "F1") is built as the charged-current (bare)
  /// form factor; p_model must then carry BOSON 24. VA_N_NPi arranges both.
  /// Parameters: GA (1.2754), fpi (0.093), gpiNN (GT value).
  /// To add a model: add an enum entry, a case in ValidateModel(), a
  /// self-contained function, and a case in FF_Assigner().
  //////////////////////////////////////////////////////////////////////////
  class FF_N_NPi : public FormFactor_Base {
  public:
    enum class Model : int {
      // none sets all to 0, pointlike sets FpiNN=g_piNN, FPiPole=1
      None = 0, Point = 1,
      // Standalone 10-19
      PionDipole = 10,
      Rho = 11,
      F1V = 12,
      // PCAC from FF_N_N 20-29 (same numbering as axial in FF_N_N)
      AxialDipole = 20,
      Meyer = 21,
      BBBA07 = 22,
      MINERvA = 23,
      Lattice = 24,
      DjukanovicLattice = 25
    };
  private:
    enum class FFtype { FpiNN, FPiPole };
    FFtype   m_fftype;
    double   m_mN, m_gA, m_fpi, m_gpiNN, m_mrho2;
    FF_N_N * p_NN;  // internal FF_N_N: FA for PCAC, F1 for F1V

    Model GetModel() const { return ModelAs<Model>(); }
    void ValidateModel() const;

    // standalone models
    double PionDipole(const double & x) const;
    double Rho(const double & x) const;
  public:
    FF_N_NPi(const FF_Parameters & params);
    ~FF_N_NPi();
    double  FF_Assigner(const double & x) const;
    Complex operator()(const ATOOLS::Vec4D_Vector& moms) override;
  };
}

#endif
