#ifndef SHERPA_Tools_Output_HepMC3_Short_H
#define SHERPA_Tools_Output_HepMC3_Short_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC3_Interface.H"
#include <iostream>
#include <memory>

namespace HepMC3 { 
  class GenRunInfo; 
  class Writer;
}
namespace HepMC = HepMC3;

namespace SHERPA {

  class Output_HepMC3_Short : public Output_Base {
    HepMC3_Interface m_hepmc3;
    HepMC::Writer *p_writer = nullptr;
    std::shared_ptr<HepMC::GenRunInfo> m_run_info;
    std::string m_basename;
    int m_iotype;
    double m_xs, m_err;

  public:
    Output_HepMC3_Short(const Output_Arguments &args);
    ~Output_HepMC3_Short();
    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();
  };

}

#endif
