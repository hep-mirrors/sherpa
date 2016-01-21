#ifndef SHRIMPS_Tools_MinBias_Parameters_H
#define SHRIMPS_Tools_MinBias_Parameters_H

#include "SHRiMPS/Tools/Parameter_Structures.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Form_Factor;
  class Omega_ik;
  
  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class MinBias_Parameters {
  private:
    run_mode::code        m_runmode;
    size_t                m_NGWstates;
    double                m_originalY, m_bmax, m_accu;
    FormFactor_Parameters m_ff_params;
    Eikonal_Parameters    m_eik_params;
    Run_Parameters        m_run_params;
    Ladder_Parameters     m_ladder_params; 
    ShowerLink_Parameters m_showerlink_params; 

    std::list<Form_Factor *> * p_ffs;
    std::list<Omega_ik *>    * p_eikonals;

    
    void FillRunParameters(ATOOLS::Data_Reader * dr);
    void FillFormFactorParameters(ATOOLS::Data_Reader * dr);
    void FillEikonalParameters(ATOOLS::Data_Reader * dr);
    void FillLadderParameters(ATOOLS::Data_Reader * dr);
    void FillShowerLinkParameters(ATOOLS::Data_Reader * dr);
  public:
    MinBias_Parameters();
    ~MinBias_Parameters();

    void Init(ATOOLS::Data_Reader * dr);
    void Reset();
    void AddFormFactor(Form_Factor * ff) { p_ffs->push_back(ff); }  
    void AddEikonal(Omega_ik * eikonal)  { p_eikonals->push_back(eikonal); }  
    
    const run_mode::code & RunMode()   const { return m_runmode; }
    const size_t         & NGWStates() const { return m_NGWstates; }
    
    const Run_Parameters & GetRunParameters() const { 
      return m_run_params; 
    }
    const FormFactor_Parameters & GetFFParameters() const { 
      return m_ff_params; 
    }
    const Eikonal_Parameters & GetEikonalParameters() const { 
      return m_eik_params; 
    }
    const Ladder_Parameters & GetLadderParameters() const { 
      return m_ladder_params; 
    }
    const ShowerLink_Parameters & GetShowerLinkParameters() const { 
      return m_showerlink_params; 
    }
    std::list<Form_Factor *> * GetFormFactors() const {
      return p_ffs;
    }
    std::list<Omega_ik *> * GetEikonals() const {
      return p_eikonals;
    }
  };
  
  extern MinBias_Parameters MBpars;
}

#endif
