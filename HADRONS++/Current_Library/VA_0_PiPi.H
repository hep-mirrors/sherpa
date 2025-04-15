#ifndef HADRONS_Current_Library_VA_0_PiPi_H
#define HADRONS_Current_Library_VA_0_PiPi_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/Resonance.H"
#include "HADRONS++/Main/Tools.H"
#include <list>


namespace HADRONS {  
  class VA_0_PiPi : public Current_Base {
  private:
    enum class PSmode {
      pipi       = 1,
      KK         = 2,
      etapi      = 3,
      etaprimepi = 4,
      Kpi        = 11,
      Keta       = 12,
      Ketaprime  = 13,
      unknown    = 99
    };
    
    enum class ffmodel {
      none  = 0,
      KS    = 1,
      RChT  = 2
    };

  public:
    static ResDefaults s_KSdefaults, s_RChTdefaults;
    
  private:
    resonance_type m_restype;
    PSmode         m_PSmode;
    ffmodel        m_ffmodel;
    double         m_norm, m_deltaM2;
    double         m_fpi, m_fa0, m_GV, m_GS, m_RChTpref;
    double         m_g_a0_eta_pi, m_g_a0_etaprime_pi;
    double         m_eps_etapi, m_eps_etaprimepi;
    double         m_m2_pi, m_m2_K, m_m2_eta, m_mu2;

    std::list<Resonance_Base * > m_vectors, m_scalars;

    void FixMode();
    void FixNorm(struct GeneralModel & model);
    bool SelectResonances(std::list<kf_code> & tags,
			  const bool & isScalar=false);
    Resonance_Base * InitResonance(struct GeneralModel & model,const kf_code & tag);
    bool FillDefaults(const kf_code & tag, std::vector<double> & defs);
    bool KSDefaults(const kf_code & tag, std::vector<double> & defs);
    bool RChTDefaults(const kf_code & tag, std::vector<double> & defs);
    
    Complex VectorFF(const double & s);
    Complex VectorKS(const double & s);
    Complex VectorRChT(const double & s);
    Complex VectorRChT_pipi(const double & s);
    Complex VectorRChT_Kpi(const double & s);
    Complex ScalarFF(const double & s);
    Complex ScalarKS(const double & s);
    Complex ScalarRChT(const double & s);
    Complex ScalarRChT_Kpi(const double & s);

    // Exponent for RChT pipi form factors.
    Complex Loop(const double & m2,const double & s,const double & mu2);
    // Exponent for RChT Kpi form factors.
    Complex Htilde(const double & s,const double & MP2,const double & MQ2);
    // Functions below from Gasser & Leutwyler, NPB250 (1985) 465:
    // JBar, JBarBar: Eq (8.8), L, MR, k: Eq (8.9)
    Complex JBar(const double & s,const double & MP2,const double & MQ2);
    Complex JBarBar(const double & s,const double & MP2,const double & MQ2);
    Complex L(const double & s,const double & MP2,const double & MQ2);
    Complex Mr(const double & s,const double & MP2,const double & MQ2);
    Complex k(const double & MP2,const double & MQ2);
  public:
    VA_0_PiPi(const ATOOLS::Flavour_Vector& flavs,
	      const std::vector<int>& indices,
	      const std::string& name);
    ~VA_0_PiPi();
    
    void SetModelParameters(struct GeneralModel model);
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}

#endif
