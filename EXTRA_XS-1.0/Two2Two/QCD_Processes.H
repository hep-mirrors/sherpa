#ifndef QCD_Process_H
#define QCD_Process_H

#include "XS_Group.H"

namespace AMEGIC { class  Process_Group; }

namespace EXTRAXS {

  class QCD_Processes: public XS_Group {
    double aS;
    double ycut;
    int    jetfinder;

    void Fill4qmodes();
    void Fill2q2gmodes();
    void Fill4gmodes();

    void CreateSelector();
  public:
    QCD_Processes();
    
    double Scale(AMATOOLS::Vec4D *);
    double KFactor(double);

    void Initialize(APHYTOOLS::Flavour *) { }
    void Initialize(ISR::ISR_Handler *, BEAM::Beam_Spectra_Handler *,
		    APHYTOOLS::Selector_Data *, AMEGIC::Process_Group *);
  };
}

#endif
