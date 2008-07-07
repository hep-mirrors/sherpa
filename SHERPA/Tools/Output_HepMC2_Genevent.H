#ifndef Output_HepMC2_Genevent_H
#define Output_HepMC2_Genevent_H

#include "Output_Base.H"
#include "HepMC2_Interface.H"

namespace HepMC { class IO_GenEvent; }

namespace SHERPA {

  class Output_HepMC2_Genevent : public Output_Base {
    HepMC2_Interface m_hepmc2;
    HepMC::IO_GenEvent* p_iogenevent;
  public:
    Output_HepMC2_Genevent(std::string basename,std::string ext,int precision);
    ~Output_HepMC2_Genevent();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile(std::string number);
  };

}

#endif
