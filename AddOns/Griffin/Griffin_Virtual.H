#ifndef Griffin_Virtual_H
#define Griffin_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "AddOns/Griffin/Griffin_Interface.H"

namespace Griffin {

  class Griffin_Virtual : public PHASIC::Virtual_ME2_Base {

  public:

    Griffin_Virtual(const PHASIC::Process_Info& pi,
                     const ATOOLS::Flavour_Vector& flavs,
                     const GriffinProcess& proc);

    ~Griffin_Virtual() override = default;

    void Calc(const ATOOLS::Vec4D_Vector& momenta) override;

  private:
    //! This process's flavours and normalisation, fixed at construction by
    //! MakeProcess() -- see the rationale on GriffinProcess and
    //! Griffin_Interface for why this is per-instance rather than shared.
    GriffinProcess m_process;
  };

}

#endif
