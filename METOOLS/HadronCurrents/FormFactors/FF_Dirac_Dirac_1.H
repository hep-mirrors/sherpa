
#ifndef METOOLS_HadronCurrents_FormFactors_FF_Dirac_Dirac_1_H
#define METOOLS_HadronCurrents_FormFactors_FF_Dirac_Dirac_1_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {

  enum class FF_Dirac_Dirac_1_mode {
    pp      = 1,
    nn      = 2,
    unknown = 999
  };

  class FF_Dirac_Dirac_1_Base : public FormFactor_Base {
  public:
    enum class Model : int {
      None      = 0,
      pp_Dipole = 400,
      pp_Kelly  = 401,
      pp_BBBA05 = 402,
      pp_YAHL18 = 403,
      nn_Dipole = 404,
      nn_Kelly  = 405,
      nn_BBBA05 = 406,
      nn_YAHL18 = 407,
      pp_Point  = 408,
      nn_Point  = 409
    };

  protected:
    Model GetModel() const {
      return ModelAs<Model>();
    }

    void ValidateModel() const {
      switch (GetModel()) {
        case Model::None:
        case Model::pp_Dipole:
        case Model::pp_Kelly:
        case Model::pp_BBBA05:
        case Model::pp_YAHL18:
        case Model::nn_Dipole:
        case Model::nn_Kelly:
        case Model::nn_BBBA05:
        case Model::nn_YAHL18:
        case Model::pp_Point:
        case Model::nn_Point:
          return;
      }

      THROW(fatal_error,
            "Unknown FORM_FACTOR=" + std::to_string(m_ffmodel) +
            " for FF_Dirac_Dirac_1_Base");
    }

    FF_Dirac_Dirac_1_mode m_mode;

    void FixMode();

    void Construct(const FF_Parameters& params);
    void Construct_pp(const FF_Parameters& params);
    void Construct_nn(const FF_Parameters& params);

    Complex F1_pp_Dipole(const double& Q2);
    Complex F1_pp_Kelly(const double& Q2);
    Complex F1_pp_BBBA05(const double& Q2);
    Complex F1_pp_YAHL18(const double& Q2);

    Complex F1_nn_Dipole(const double& Q2);
    Complex F1_nn_Kelly(const double& Q2);
    Complex F1_nn_BBBA05(const double& Q2);
    Complex F1_nn_YAHL18(const double& Q2);

    double Mp;
    double mu_p;
    double Mn;
    double mu_n;
    double Lambda2;

  public:
    explicit FF_Dirac_Dirac_1_Base(const FF_Parameters& params)
      : FormFactor_Base(params),
        m_mode(FF_Dirac_Dirac_1_mode::unknown)
    {
      ValidateModel();
      FixMode();
      Construct(params);
      CreateF1Histogram();  
    }

    virtual ~FF_Dirac_Dirac_1_Base() = default;

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);

    void CreateF1Histogram();
  };

}

#endif

