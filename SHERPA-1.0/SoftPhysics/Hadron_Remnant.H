#ifndef Hadron_Remnant_H
#define Hadron_Remnant_H

#include "QCD_Remnant_Base.H"

namespace SHERPA {

  class Hadron_Remnant: public QCD_Remnant_Base {
  private:

    ATOOLS::Vec4D m_hardpt;

    size_t m_last;

    void GetConstituents(ATOOLS::Flavour flav);

    double GetXPDF(ATOOLS::Flavour flavour,double scale);
    void   DiceKinematics();

    bool AttachLastRemnants();

  public:
    
    // constructor
    Hadron_Remnant(PDF::ISR_Handler *isrhandler,
		   const unsigned int beam,const double scale);

    // member functions
    double MinimalEnergy(const ATOOLS::Flavour &flavour); 
    
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);

    void UnDo();

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
