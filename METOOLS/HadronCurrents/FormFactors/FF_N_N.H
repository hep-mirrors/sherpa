#ifndef METOOLS_Current_Library_FF_N_N_H
#define METOOLS_Current_Library_FF_N_N_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  // Boson x nucleon line, fixed from the legs and BOSON in FixMode().
  enum class FF_N_N_vertex {
    gamma_p = 1, gamma_n = 2,
    Z_p     = 11, Z_n    = 12,
    W_pn    = 21,
    unknown = 999
  };

  //////////////////////////////////////////////////////////////////////////
  /// Nucleon form factors for VA_N_N, already combined for the boson:
  ///   F1, F2 : photon F^N, Z T3 F^V - 2 s_W^2 F^N, W F^V = F^p - F^n
  ///   FA     : photon 0, Z T3 F_A, W F_A
  ///   FP     : PCAC, 4 m_N^2/(Q^2+m_pi^2) FA
  /// FORM_FACTOR: 0 zero, 1 point (F1p=1, F_A=g_A),
  ///   vector 10-19 : 10 Kelly
  ///   axial  20-29 : 20 Meyer deuterium z-exp., 21 lattice z-exp.
  /// To add a model: add an enum entry, a case in ValidateModel(), a
  /// self-contained function, and a case in Vector() or Axial().
  //////////////////////////////////////////////////////////////////////////
  class FF_N_N : public FormFactor_Base {
  public:
    enum class Model : int {
      // none sets all to 0, pointlike sets F1p=1, F_A=g_A, all else 0
      None = 0, Point = 1,
      // vector 10-19
      Dipole = 10,
      Kelly = 11,
      BBBA05 = 12,
      Ye2018 = 13,
      // axial 20-29
      AxialDipole = 20,
      Meyer = 21,
      BBBA07 = 22,
      MINERvA = 23,
      Lattice = 24,
      DjukanovicLattice = 25
    };
  protected:
    Model GetModel() const { return ModelAs<Model>(); }
    void ValidateModel() const {
      const bool axial(m_fftype==FFtype::FA || m_fftype==FFtype::FP);
      switch (GetModel()) {
        case Model::None: case Model::Point:
          return;
        case Model::Dipole:
        case Model::Kelly:
        case Model::BBBA05:
        case Model::Ye2018:

        case Model::AxialDipole:
        case Model::Meyer:
        case Model::BBBA07:
        case Model::MINERvA:
        case Model::Lattice:
        case Model::DjukanovicLattice:
          if (m_ffmodel/10==(axial ? 2 : 1)) return;
      }
      THROW(fatal_error, "Unknown FORM_FACTOR=" + std::to_string(m_ffmodel) +
            " for FF_N_N '" + m_name + "'");
    }
  private:
    enum class FFtype { F1, F2, FA, FP };
    FFtype       m_fftype;
    FF_N_N_vertex m_vertex;
    double      m_mN, m_mpi, m_s2, m_T3, m_gA;

    void FixMode(const GeneralModel & md);

    // vector models: each returns F1, F2 of one nucleon
    void   Dipole_p(const double & Q2,double & F1,double & F2) const;
    void   Dipole_n(const double & Q2,double & F1,double & F2) const;
    void   Kelly_p(const double & Q2,double & F1,double & F2) const;
    void   Kelly_n(const double & Q2,double & F1,double & F2) const;
    void   BBBA05_p(const double & Q2,double & F1,double & F2) const;
    void   BBBA05_n(const double & Q2,double & F1,double & F2) const;
    void   Ye2018_p(const double & Q2,double & F1,double & F2) const;
    void   Ye2018_n(const double & Q2,double & F1,double & F2) const;
    // axial models: each returns the normalised F_A
    double AxialDipole(const double & Q2) const;
    double Meyer(const double & Q2) const;
    double BBBA07(const double & Q2) const;
    double MINERvA(const double & Q2) const;
    double Lattice(const double & Q2) const;
    double DjukanovicLattice(const double & Q2) const;

    void   Vector(const double & Q2,double & F1p,double & F2p,
                  double & F1n,double & F2n) const;
    double Axial(const double & Q2) const;
  public:
    FF_N_N(const FF_Parameters & params);
    ~FF_N_N() {}
    double  FF_Assigner(const double & Q2) const;
    void CreateF1F2Histogram() const;
    void CreateAxialRatioHistogram() const;
    Complex operator()(const ATOOLS::Vec4D_Vector& moms) override;
  };
}

#endif
