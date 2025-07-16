#ifndef HADRONS_ME_Library_Two_Body_MEs_H
#define HADRONS_ME_Library_Two_Body_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class Baryon_Diquark_Quark : public HD_ME_Base {
    double m_Vxx_decay, m_Vxx_production;
    double m_GF;
  public:
    Baryon_Diquark_Quark(const ATOOLS::Flavour_Vector& flavs,
                         const std::vector<int>& decayindices,
                         const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
    bool SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q,int n_g,
                      bool m_anti);
  };

} // end of namespace


#endif
