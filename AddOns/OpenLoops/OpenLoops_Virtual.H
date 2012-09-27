#ifndef OpenLoops_Virtual_H
#define OpenLoops_Virtual_H

#include "OpenLoops_Interface.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include <iostream>
#include <fstream>

namespace OpenLoops {

  class OpenLoops_Virtual : public PHASIC::Virtual_ME2_Base {

    Amp2Func m_amp2;
    PermutationFunc m_permutationfunc;
    std::vector<int> m_permutation;
    static OpenLoops_Interface* s_interface;
    int m_oew, m_oqcd;
    int m_write_points;
    std::ofstream m_points_file;

    bool     m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

  public:

    OpenLoops_Virtual(const PHASIC::Process_Info& pi,
                      const ATOOLS::Flavour_Vector& flavs,
                      Amp2Func amp2,
                      PermutationFunc permutationfunc,
                      std::vector<int> permutation,
                      std::string functag);

    ~OpenLoops_Virtual();
  
    inline bool SetColours(const ATOOLS::Vec4D_Vector& momenta) {
      return true;
    }

    inline double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  
    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    inline static void SetInterface(OpenLoops_Interface* interface)
    {
      s_interface=interface;
    }

  };

}

#endif
