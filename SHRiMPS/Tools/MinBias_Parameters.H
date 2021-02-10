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
  class Cross_Sections;

  struct XSecs_Container {
    double xs_tot, xs_in, xs_el, xs_SD[2], xs_DD;
    XSecs_Container(): xs_tot(0.), xs_in(0.), xs_el(0.), xs_DD(0.)
    { xs_SD[0] = xs_SD[1] = 0.; }
  };
  
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

    std::vector<Form_Factor *>            m_ffs;
    std::vector<std::vector<Omega_ik *> > m_eikonals;
    XSecs_Container                       m_xsecs;
    
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
    void ResetEikonals(const size_t size);
    void AddFormFactor(Form_Factor * ff) { m_ffs.push_back(ff); }  
    void AddEikonal(const size_t & i,const size_t & j,
		    Omega_ik * eikonal)  { m_eikonals[i][j] = eikonal; }  
    void SetXSecs(Cross_Sections * xsecs);
    void UpdateForNewEnergy(const double & energy);
    
    const run_mode::code & RunMode()   const { return m_runmode; }
    const size_t         & NGWStates() const { return m_NGWstates; }
    
    const Run_Parameters        & GetRunParameters()        const { return m_run_params; }
    const FormFactor_Parameters & GetFFParameters()         const { return m_ff_params; }
    const Eikonal_Parameters    & GetEikonalParameters()    const { return m_eik_params; }
    const Ladder_Parameters     & GetLadderParameters()     const { return m_ladder_params; }
    const ShowerLink_Parameters & GetShowerLinkParameters() const { return m_showerlink_params; }
    std::vector<Form_Factor *>            * GetFormFactors() { return &m_ffs; }
    std::vector<std::vector<Omega_ik *> > * GetEikonals()    { return &m_eikonals; }
    XSecs_Container                       * GetXSecs()       { return &m_xsecs; } 
  };
  
  extern MinBias_Parameters MBpars;
}

#endif
