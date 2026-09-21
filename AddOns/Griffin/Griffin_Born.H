#ifndef AddOns_Griffin_Griffin_Born_H
#define AddOns_Griffin_Griffin_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/Griffin/Griffin_Interface.H"

namespace MODEL {
    class Model_Base;
}

namespace Griffin {

    class Griffin_Born: public PHASIC::Tree_ME2_Base {

    public:

    Griffin_Born(const PHASIC::External_ME_Args& args, const GriffinProcess& proc);

    ~Griffin_Born() override = default;

    double Calc(const ATOOLS::Vec4D_Vector& momenta) override;

    private:
    //! This process's flavours and normalisation, fixed at construction by
    //! MakeProcess() -- see the rationale on GriffinProcess and
    //! Griffin_Interface for why this is per-instance rather than shared.
    GriffinProcess m_process;
    };

}

#endif
