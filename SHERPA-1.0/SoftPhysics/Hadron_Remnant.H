#ifndef Hadron_Remnant_H
#define Hadron_Remnant_H

#include "QCD_Remnant_Base.H"

namespace SHERPA {

  class Hadron_Remnant: public QCD_Remnant_Base {
  private:

    ATOOLS::Vec4D m_hardpt;

    void GetConstituents(ATOOLS::Flavour flav);

    double GetXPDF(ATOOLS::Flavour flavour,double scale);
    void   DiceKinematics();

    bool TreatFirstQuark(ATOOLS::Particle *cur);
    bool TreatFirstGluon(ATOOLS::Particle *cur);
    bool TreatQuark(ATOOLS::Particle *cur);
    bool TreatGluon(ATOOLS::Particle *cur);

  public:
    
    // constructor
    Hadron_Remnant(PDF::ISR_Handler *isrhandler,
		   const unsigned int beam,const double scale);

    // member functions
    double MinimalEnergy(const ATOOLS::Flavour &flavour); 
    
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
