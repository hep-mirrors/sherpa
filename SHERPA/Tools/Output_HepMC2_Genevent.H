#ifndef SHERPA_Tools_Output_HepMC2_Genevent_H
#define SHERPA_Tools_Output_HepMC2_Genevent_H

#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC2_Interface.H"

namespace HepMC { class IO_GenEvent; class GenCrossSection; }

namespace SHERPA {

  class Output_HepMC2_Genevent : public Output_Base {
    HepMC2_Interface m_hepmc2;
    HepMC::IO_GenEvent* p_iogenevent;
    HepMC::GenCrossSection* p_xs;
  public:
    Output_HepMC2_Genevent(std::string basename,std::string ext,int precision);
    ~Output_HepMC2_Genevent();
    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile(std::string number);
  };

}

#endif
