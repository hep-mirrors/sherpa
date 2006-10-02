#ifndef Four_Body_MEs_H
#define Four_Body_MEs_H

#include "HD_ME_Base.H"

namespace HADRONS {

  class Higgs_TauTau_2Pi2Neutrino : public HD_ME_Base {
      bool m_pseudo;
    public:
      Higgs_TauTau_2Pi2Neutrino( int nout, ATOOLS::Flavour *fl );
      void   operator()(
                        const ATOOLS::Vec4D  * _p,
                        ATOOLS::Amplitude_Tensor     * amps,
                        int                    k0_n
                       );
      void   SetModelParameters( GeneralModel _md );
  };

}


#endif
