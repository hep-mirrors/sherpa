#ifndef METOOLS_Current_Library_FF_N_Nprime_H
#define METOOLS_Current_Library_FF_N_Nprime_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  // Boson x nucleon line, fixed from the legs and BOSON in FixMode().
  enum class FF_N_Nprime_vertex {
    gamma_p = 1, gamma_n = 2,
    unknown = 999
  };

//////////////////////////////////////////////////////////////////////////
// Electromagnetic N -> N* transition form factors.
//
// F1prime, F2prime: transition Dirac and Pauli form factors.
// Proton and neutron transitions are selected using m_vertex.
//
// FORM_FACTOR: 0 zero, 1 point.
// To add a model: add an enum entry, a case in ValidateModel(),
// and implement the model in Vector().
//////////////////////////////////////////////////////////////////////////

  class FF_N_Nprime : public FormFactor_Base {
  public:
    enum class Model : int {
      // none sets all to 0, pointlike sets F1p=1, F_A=g_A, all else 0
      None = 0, Point = 1,
      // N1440 10-19
      MAID2009_1440 = 10,

      // N1535 20-29
      MAID2009_1535 = 20


    };
  protected:
    Model GetModel() const { return ModelAs<Model>(); }
    void ValidateModel() const {
      switch (GetModel()) {
        case Model::None:
        case Model::Point:
        case Model::MAID2009_1440:
        case Model::MAID2009_1535:
          return;
      }
      THROW(fatal_error, "Unknown FORM_FACTOR=" + std::to_string(m_ffmodel) +
            " for FF_N_Nprime '" + m_name + "'");
    }
  private:
    enum class FFtype { F1prime, F2prime };
    FFtype       m_fftype;
    FF_N_Nprime_vertex m_vertex;
    double      m_N, m_R;

    void FixMode(const GeneralModel & md);

    void   MAID2009_1440_p(const double & Q2,double & F1,double & F2) const;
    void   MAID2009_1440_n(const double & Q2,double & F1,double & F2) const;
    void   MAID2009_1535_p(const double & Q2,double & F1,double & F2) const;
    void   MAID2009_1535_n(const double & Q2,double & F1,double & F2) const;

    void   Vector(const double & Q2,double & F1primep,double & F2primep,
                  double & F1primen,double & F2primen) const;
  public:
    FF_N_Nprime(const FF_Parameters & params);
    ~FF_N_Nprime() {}
    double  FF_Assigner(const double & Q2) const;
    void CreateF1F2Histogram() const;
    Complex operator()(const ATOOLS::Vec4D_Vector& moms) override;
  };
}

#endif
