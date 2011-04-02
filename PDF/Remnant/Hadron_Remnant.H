#ifndef PDF_Remnant_Hadron_Remnant_H
#define PDF_Remnant_Hadron_Remnant_H

#include "PDF/Remnant/QCD_Remnant_Base.H"

namespace PDF {

  class Hadron_Remnant: public QCD_Remnant_Base {
  private:

    ATOOLS::Vec4D m_hardpt;

    double m_xrem;
    size_t m_last;

    bool   GenerateKinematics();

    ATOOLS::Flavour Opposite(ATOOLS::Flavour flav) const; 

    bool ValenceQuark(ATOOLS::Particle *const quark); 
    bool DecomposeHadron();

  public:
    
    // constructor
    Hadron_Remnant(PDF::ISR_Handler *isrhandler,
		   const unsigned int beam);

    // member functions
    double          MinimalEnergy(const ATOOLS::Flavour &flavour); 
    
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);

    const std::vector<ATOOLS::Flavour> &GetConstituents(ATOOLS::Flavour flav);

    double GetXPDF(ATOOLS::Flavour flavour,double scale);

  };// end of class Remnant_Base

}// end of namespace PDF

#endif
