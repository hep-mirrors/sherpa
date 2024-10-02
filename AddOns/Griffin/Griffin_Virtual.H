#ifndef Griffin_Virtual_H
#define Griffin_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Base.H"

namespace Griffin {

  class Griffin_Virtual : public PHASIC::Virtual_ME2_Base {

  public:
 
    Griffin_Virtual(const PHASIC::Process_Info& pi,
                   const ATOOLS::Flavour_Vector& flavs);

    ~Griffin_Virtual() {}


    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };

}

#endif
