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
                     const GriffinProcess& proc,
                     const bool higherorder = false);

    ~Griffin_Virtual() override = default;

    void Calc(const ATOOLS::Vec4D_Vector& momenta) override;
    //! GriffinHO only: NNLO - NLO per helicity at the last Calc point
    //! (Griffin_Interface::HigherOrderHelicity)
    bool HelicityAmplitudes(std::vector<Complex> &a0, std::vector<Complex> &a1,
                            std::vector<std::vector<int> > &hel) override;

  private:
    //! This process's flavours and normalisation, fixed at construction by
    //! MakeProcess() -- see the rationale on GriffinProcess and
    //! Griffin_Interface for why this is per-instance rather than shared.
    GriffinProcess m_process;
    //! Loop generator "GriffinHO": return NNLO - NLO (YFS: EW_HIGHER_ORDER)
    //! instead of GRIFFIN: Order - LO.
    bool m_higherorder;
    ATOOLS::Vec4D_Vector m_lastmom;
    double m_norm = 1./M_PI;
  };

}

#endif
