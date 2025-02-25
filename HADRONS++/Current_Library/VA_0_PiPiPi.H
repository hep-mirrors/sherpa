#ifndef HADRONS_Current_Library_VA_0_PiPiPi_H
#define HADRONS_Current_Library_VA_0_PiPiPi_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/Resonance.H"
#include "HADRONS++/Main/Tools.H"
#include <list>


namespace HADRONS {  
  class VA_0_PiPiPi : public Current_Base {
  private:
    enum class PSmode {
      pipipi_pmp       = 1,
      pipipi_p00       = 2,
      unknown          = 99
    };
    
    enum class ffmodel {
      none = 0,
      KS   = 1
    };
        
    resonance_type m_restype;
    PSmode         m_PSmode;
    ffmodel        m_ffmodel;
    double  m_global, m_deltaM2;
    double  m_fpi, m_m2_pi, m_mu2;

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
