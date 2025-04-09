#ifndef HADRONS_Current_Library_VA_0_PiPiPi_H
#define HADRONS_Current_Library_VA_0_PiPiPi_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/Resonance.H"
#include "HADRONS++/Main/Tools.H"
#include <list>


namespace HADRONS {  
  class VA_0_PiPiPi : public Current_Base {
  public:
    static ResDefaults s_KSdefaults;
  private:      
    enum class PSmode {
      pipipi_pmp       = 1,
      pipipi_p00       = 2,
      Kpipi_pmp        = 11,
      Kpipi_p00        = 12,
      piKpi_ppm        = 13,
      piKpi_p00        = 14,
      unknown          = 99
    };
    
    enum class ffmodel {
      none = 0,
      KS   = 1
    };

    enum class fftag {
      ffV_1 = 1,
      ffV_2 = 2,
      ffV_3 = 3,
      ffS_4 = 4,
      ffS_5 = 5
    };

    PSmode         m_PSmode;
    ffmodel        m_ffmodel;
    resonance_type m_restype;
    double         m_norm;
    double         m_fpi, m_m2_pi, m_mu2;
    size_t         m_i, m_j, m_k;
    Complex        m_c1, m_c2, m_c4, m_c5;
    
    std::map<fftag,Resonance_Tree * >           m_treemap;
    std::map<ATOOLS::Flavour,Resonance_Base * > m_initialised;

    void    FixMode();
    void    FixMappings();
    void    FixNorm(struct GeneralModel & model);
    bool    SelectResonances(struct GeneralModel & model);
    void    MakeTree(struct GeneralModel & model,fftag ff,
		     const ATOOLS::Flavour & start,
		     const std::vector<ATOOLS::Flavour> & ijk,
		     const std::list<ATOOLS::Flavour> & rlist,
		     const std::vector<ATOOLS::Flavour> & ij);
    Resonance_Base * MakeResonance(struct GeneralModel & model,
				   const ATOOLS::Flavour & flav,
				   const std::vector<ATOOLS::Flavour> & ij);
    bool    FillDefaults(const kf_code & tag, std::vector<double> & defs);
    bool    KSDefaults(const kf_code & tag, std::vector<double> & defs);
    Complex FF_V1(const double & sijk,const double & sij);
    Complex FF_V2(const double & sijk,const double & skj);
    Complex FF_A(const double & s);
    Complex FF_P(const double & s);
  public:
    VA_0_PiPiPi(const ATOOLS::Flavour_Vector& flavs,
		const std::vector<int>& indices,
		const std::string& name);
    ~VA_0_PiPiPi();
    
    void SetModelParameters(struct GeneralModel model);
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}

#endif
