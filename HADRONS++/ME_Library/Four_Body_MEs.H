#ifndef HADRONS_ME_Library_Four_Body_MEs_H
#define HADRONS_ME_Library_Four_Body_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class QQ_QQQQ_Spectator : public HD_ME_Base {
    ATOOLS::Flavour m_decayer;
    double m_Vxx_decay, m_Vxx_production;
    double m_GF;
    double m_cR_decay, m_cL_decay, m_cR_production, m_cL_production;
    double m_colourflip_ratio;
  public:
    QQ_QQQQ_Spectator(const ATOOLS::Flavour_Vector& flavs,
		      const std::vector<int>& decayindices,
		      const std::string& name) :
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
    bool SetColorFlow(std::vector<ATOOLS::Particle*> outparts,
			int n_q, int n_g, bool anti);
  };
  
  class Baryon_Diquark_3Quarks : public HD_ME_Base {
    double m_Vxx_decay, m_Vxx_production;
    double m_GF;
  public:
    Baryon_Diquark_3Quarks(const ATOOLS::Flavour_Vector& flavs,
                           const std::vector<int>& decayindices,
                           const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void   SetModelParameters( GeneralModel _md );
    bool   SetColorFlow(std::vector<ATOOLS::Particle*> outparts,
                        int n_q, int n_g, bool anti);
  };
  
  class B_tautau_pinupinu : public HD_ME_Base {
  public:
    B_tautau_pinupinu(const ATOOLS::Flavour_Vector& flavs,
                      const std::vector<int>& decayindices,
                      const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };
}


#endif
