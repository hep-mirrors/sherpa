#ifndef MODEL_DarkQCD_Model_H
#define MODEL_DarkQCD_Model_H

#include "MODEL/SM/Model.H"
#include "MODEL/Main/Running_AlphaS.H"

namespace HADRONS {
  class Hadron_Decay_Table;
  class Mixing_Handler;
  struct GeneralModel;
};

namespace ATOOLS {
  class Mass_Selector;
};

namespace MODEL {

  class One_Running_AlphaDark : public One_Running_AlphaS {
  private:
    double m_Nc;
    void SetNC(double NC) {
      m_Nc = NC;
      m_TR = 1./2.;
      m_CF = m_TR*(m_Nc*m_Nc-1.)/(m_Nc);
      m_CA = 2*m_TR*m_Nc;
    }
  public:

    One_Running_AlphaDark(const double as_MZ = 0.0,const double m2_MZ = 0.0,
                          const int order = 0, const int thmode = 1) :
      One_Running_AlphaS(NULL,as_MZ,m2_MZ,order,thmode)
    {
      auto s = ATOOLS::Settings::GetMainSettings()["DARK_SECTOR"];
      SetNC(s["NC"].SetDefault(3).Get<double>());
      /// TODO: Set thresholds.
    }
    
  };

  class Running_AlphaDark : public Running_AlphaS {
  public:
    Running_AlphaDark(const double as_MZ,const double m2_MZ,
                      const int order, const int thmode,
                      const PDF::ISR_Handler_Map& isr) : Running_AlphaS(as_MZ,m2_MZ,
                                                                        order,thmode, isr) {
      m_type = "Running_Coupling";
      m_name = "Alpha_Dark";
      m_alphas.clear();
      for (PDF::ISR_Handler_Map::const_iterator it=isr.begin(); it!=isr.end(); ++it) {
        m_alphas.insert(std::make_pair(it->first, new One_Running_AlphaDark(as_MZ, m2_MZ, order, thmode)));
      }
      SetActiveAs(PDF::isr::hard_process); 
    }
    
  };

  
  class DarkQCD: public Standard_Model {
  private:
    void ParticleZprimeInit();
    void InitZprimeVertices();
    void InitDarkQCDVertices();
    bool m_dec_dark_g4;

    double m_fPi, m_fEta, m_gRhoPiPi, m_gOmegaPiPi;

    const ATOOLS::Mass_Selector   * p_ms;
    HADRONS::Mixing_Handler * p_mix;
    struct SortByMass {
      bool operator()(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B) const {
	return A.Mass() > B.Mass();
      };
    };
    std::set<ATOOLS::Flavour, SortByMass> m_hadrons;
    HADRONS::Hadron_Decay_Table * DecaysForRhos(const ATOOLS::Flavour & flav);
    HADRONS::Hadron_Decay_Table * DecaysForPis(const ATOOLS::Flavour & flav);

    double GammaTwoBody(const double & Ampl2,const double & M2,
			const double & m12,const double & m22);
  public:
    DarkQCD();

    void AddHadronDecays(PHASIC::Decay_Map * dmap);
    bool ModelInit();
    void InitVertices();
    
  };

}

#endif
