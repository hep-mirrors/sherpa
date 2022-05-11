#ifndef PHASIC_Scales_Pilot_Scale_Setter_H
#define PHASIC_Scales_Pilot_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/METS_Scale_Setter.H"
#include "PHASIC++/Scales/Variable_Scale_Setter.H"

namespace PHASIC {

  class Pilot_Scale_Setter : public Scale_Setter_Base {

  public:

    Pilot_Scale_Setter(const Scale_Setter_Arguments &args,
                       const int mets_mode = 2);

    double CalculateScale(const ATOOLS::Vec4D_Vector &p,const size_t mode=0) override;
    PDF::CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const override;
    void SetCaller(Process_Base *const proc) override;
    double Scale(const ATOOLS::stp::id type,const int mode=0) const override;
    const std::vector<double> &Scales() const override;
    const std::vector<double> &FixedScales() const override;
    ATOOLS::ClusterAmplitude_Vector &Amplitudes() override;
    const ATOOLS::Vec4D_Vector &Momenta() const override;
    void PreCalc(const ATOOLS::Vec4D_Vector &p, const size_t &mode) override;
    double Calculate(const ATOOLS::Vec4D_Vector &p,
                     const size_t &mode) override;

  private:
    METS_Scale_Setter mets;
    Variable_Scale_Setter var;

  };

}

#endif
