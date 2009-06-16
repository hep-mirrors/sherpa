#ifndef HADRONS_ME_Library_Two_Body_MEs_H
#define HADRONS_ME_Library_Two_Body_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class Baryon_Diquark_Quark : public HD_ME_Base {
    double m_Vxx_decay, m_Vxx_production;
    double m_GF;
  public:
    Baryon_Diquark_Quark(ATOOLS::Flavour* flavs, int n, int* decayindices,
                         std::string name):
      HD_ME_Base(flavs,n,decayindices,name) {};
    void operator()(
		      const ATOOLS::Vec4D  * p,
		      HELICITIES::Spin_Amplitudes     * amps
		      );
    void SetModelParameters( GeneralModel _md );
    bool SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q,int n_g);
  };

} // end of namespace


#endif
