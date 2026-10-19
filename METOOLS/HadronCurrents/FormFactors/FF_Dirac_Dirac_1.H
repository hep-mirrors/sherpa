
#ifndef METOOLS_HadronCurrents_FormFactors_FF_Dirac_Dirac_1_H
#define METOOLS_HadronCurrents_FormFactors_FF_Dirac_Dirac_1_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {

  enum class FF_Dirac_Dirac_1_mode {
    pp              = 1,
    nn              = 2,
    unknown         = 999
  };


class FF_Dirac_Dirac_1_Base : public FormFactor_Base {
    protected:
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
      FF_Dirac_Dirac_1_Base(const FF_Parameters& params)
        : FormFactor_Base(params),
          m_mode(FF_Dirac_Dirac_1_mode::unknown)
      {
        FixMode();
        Construct(params);
      }
    
      virtual ~FF_Dirac_Dirac_1_Base() = default;
    
      Complex operator()(const ATOOLS::Vec4D_Vector& moms);
};

}


#endif
