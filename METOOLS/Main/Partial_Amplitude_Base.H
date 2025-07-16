#ifndef METOOLS_Main_Partial_Amplitude_Base_H
#define METOOLS_Main_Partial_Amplitude_Base_H

#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {

  class  Partial_Amplitude_Base : public Spin_Amplitudes {
  protected:
    const ATOOLS::Flavour_Vector& p_flavs;
    std::vector<int>              p_i; // external indices for flavs and moms
    std::vector<bool>             p_out; // incoming/outgoing

    void AssertSpins(int, ...);
    void AssertIn(int in);
  public:
    Partial_Amplitude_Base(const ATOOLS::Flavour_Vector& flavs,
                           const std::vector<int>& i,
                           const std::vector<bool>& out);
    virtual ~Partial_Amplitude_Base();

    static Partial_Amplitude_Base* Select(const ATOOLS::Flavour_Vector& flavs);
  };

  class Isotropic : public Partial_Amplitude_Base {
  public:
    Isotropic(const ATOOLS::Flavour_Vector& flavs,
              const std::vector<int>& i, const std::vector<bool>& out);
    ~Isotropic() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

}


#endif
