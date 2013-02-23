#ifndef AddOns_OpenLoops_OpenLoops_Born_H
#define AddOns_OpenLoops_OpenLoops_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/OpenLoops/OpenLoops_Interface.H"

namespace MODEL {
  class Model_Base;
}

namespace OpenLoops {

  class OpenLoops_Born: public PHASIC::Tree_ME2_Base {

    Amp2Func m_amp2;
    PermutationFunc m_permutationfunc;
    std::vector<int> m_permutation;
    static OpenLoops_Interface* s_interface;

    size_t m_oqcd, m_oew;
    double m_symfac;

  public:

    OpenLoops_Born(const PHASIC::Process_Info& pi,
                   const ATOOLS::Flavour_Vector& flavs,
                   Amp2Func amp2,
                   PermutationFunc permutationfunc,
                   std::vector<int> permutation,
                   std::string functag);

    ~OpenLoops_Born();

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

    inline static void SetInterface(OpenLoops_Interface* interface)
    {
      s_interface=interface;
    }
  };

}

#endif
