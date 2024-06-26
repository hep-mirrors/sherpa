#ifndef PHASIC_Main_Phase_Space_Handler_H
#define PHASIC_Main_Phase_Space_Handler_H

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Point.H"
#include "PHASIC++/Main/Phase_Space_Enhance.H"
#include "PHASIC++/Main/Phase_Space_Variation_Weights.H"

#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Phys/Weights.H"
#include <list>

namespace ATOOLS { 
  class Blob_Data_Base; 
  class Poincare;
  class Mass_Selector;
  struct Weight_Info;
}

namespace BEAM {
  class Beam_Spectra_Handler;
}

namespace PDF {
  class ISR_Handler;
}

namespace YFS {
  class YFS_Handler;
}

namespace PHASIC {
  class Multi_Channel;
  class Beam_Channels;
  class ISR_Channels;
  class Subprocess_Info;
  class Process_Info;

  class Phase_Space_Integrator;
  class Process_Base;

  class Phase_Space_Handler {
  private: 
    BEAM::Beam_Spectra_Handler * p_beamhandler;
    PDF::ISR_Handler           * p_isrhandler;
    YFS::YFS_Handler           * p_yfshandler;
    std::string                  m_name;
    int                          m_nin, m_nout, m_nvec;
    ATOOLS::Flavour_Vector       p_flavours;
    ATOOLS::Vec4D_Vector         p_lab;

    Process_Integrator          * p_process, *p_active;
    Phase_Space_Integrator      * p_integrator;
    Phase_Space_Point             m_pspoint;
    Phase_Space_Enhance           m_psenhance;

    ATOOLS::Weights_Map           m_wgtmap;
    static ATOOLS::Integration_Info * p_info;

    double m_enhanceweight, m_ISsymmetryfactor, m_psweight;

    psmode::code m_cmode;
    double    m_thkill, m_error, m_abserror;
    bool      m_initialized, m_sintegrator, m_fin_opt, m_printpspoint;
    long int  m_killedpoints;
    std::vector<std::vector<double> > m_stats;

    std::string m_enhanceObs, m_enhanceFunc;

    void RegisterDefaults() const;
    void InitParameters(const double & error);

    void ManageWeights(const double & factor);
    bool CheckStability();

    void PrintIntermediate();
    void CheckSinglePoint();
    static void TestPoint(ATOOLS::Vec4D *const p,
			  ATOOLS::Vec4D_Vector cp,ATOOLS::Flavour_Vector fl,
			  const Subprocess_Info *info,size_t &n,
                          const ATOOLS::Mass_Selector* ms);

  public:

    //constructor
    Phase_Space_Handler(Process_Integrator *,double error,
                        const std::string eobs,const std::string efunc);

    //destructor
    ~Phase_Space_Handler();

    bool CreateIntegrators();
    inline void InitCuts() { m_pspoint.InitCuts(p_process); }

    double Integrate();
    ATOOLS::Weight_Info *OneEvent(
        Process_Base *const proc,
        ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::all,
        const int mode=0);
    inline ATOOLS::Weights_Map Differential(ATOOLS::Variations_Mode varmode) {
      return Differential(p_process, varmode);
    }
    ATOOLS::Weights_Map Differential(Process_Integrator *const process,
                                     ATOOLS::Variations_Mode,
                                     const psmode::code mode=psmode::normal_call);

    inline ATOOLS::Weights_Map CalculateME(ATOOLS::Variations_Mode varmode)  {
      return p_active->Process()->Differential(p_lab, varmode);
    }
    inline double CalculatePS()  { return m_pspoint.CalculateWeight(); }
    inline double Weight(ATOOLS::Vec4D_Vector& moms) {
      m_pspoint.SetMomenta(moms);
      return m_pspoint.CalculateWeight();
    }

    void   AddPoint(const double xs);

    static void TestPoint(ATOOLS::Vec4D *const p,const Process_Info *info,
                          const ATOOLS::Mass_Selector* ms,const int mode=0);
    static void TestPoint(ATOOLS::Vec4D *const p,
			  const size_t &nin,const size_t &nout,
			  const ATOOLS::Flavour_Vector &flavs,
			  const ATOOLS::Mass_Selector* ms);


    inline static ATOOLS::Integration_Info* GetInfo() {
      if (p_info==NULL) p_info = new ATOOLS::Integration_Info();
      return p_info;
    }
    inline static void DeleteInfo() {
      if (p_info) { delete p_info; p_info=NULL; }
    }

    void   WriteOut(const std::string &path);
    bool   ReadIn(const std::string &path,const size_t exclude=0);
    bool   YFSEvent();

    inline void MPISync() {
      m_pspoint.MPISync();
      p_process->MPISync();
    }
    inline void Optimize() {
      m_pspoint.Optimize(m_error);
      p_process->ResetMax(2);
      m_psenhance.Optimize();
    }
    inline void EndOptimize() {
      m_pspoint.EndOptimize(m_error);
    }
    inline void SetEnhanceWeight(const double& enhance)
    {
      m_enhanceweight = enhance;
      m_psenhance.SetFactor(enhance);
    }
    inline void SetOSMass(const double &osmass) { m_pspoint.SetOSMass(osmass); }
    inline void SetError(const double error)    { m_error=error; }
    inline void SetAbsError(const double error) { m_abserror=error; }

    inline void SetBeamIntegrator(Multi_Channel * channels) { m_pspoint.SetBeamIntegrator(channels); }
    inline void SetISRIntegrator(Multi_Channel * channels)  { m_pspoint.SetISRIntegrator(channels); }
    inline void SetFSRIntegrator(Multi_Channel * channels)  { m_pspoint.SetFSRIntegrator(channels); }

    inline BEAM::Beam_Spectra_Handler   * GetBeamSpectra() const { return p_beamhandler; }
    inline PDF::ISR_Handler              * GetISRHandler() const { return p_isrhandler; }
    inline YFS::YFS_Handler              * GetYFSHandler() const { return p_yfshandler; }
    inline Phase_Space_Integrator           * Integrator() const { return p_integrator; }
    inline const std::vector<std::vector<double> > Stats() const { return m_stats; }
    inline ATOOLS::Vec4D_Vector &                Momenta()       { return p_lab; }
    inline ATOOLS::Flavour_Vector &                Flavs()       { return p_flavours; }
    inline Cut_Data                               * Cuts()       { return m_pspoint.Cuts(); }
    inline Process_Integrator *                  Process() const { return p_process; }
    inline Process_Integrator *            Active()  const { return p_active;  }
    inline const double &                 PSWeight() const { return m_psweight; }
    inline const double &            EnhanceWeight() const { return m_enhanceweight;  }
    inline const double &                    Error() const { return m_error; }
    inline const double &                 AbsError() const { return m_abserror; }
    inline Multi_Channel *          BeamIntegrator() const { return m_pspoint.BeamIntegrator(); }
    inline Multi_Channel *           ISRIntegrator() const { return m_pspoint.ISRIntegrator(); }
    inline Multi_Channel *           FSRIntegrator() const { return m_pspoint.FSRIntegrator(); }
    inline void AddStats(const std::vector<double> &stats) {
      std::vector<double> nstats(1,m_stats.size()+1);
      nstats.insert(nstats.end(),stats.begin(),stats.end());
      m_stats.push_back(nstats); 
    }
  };// end of class Phase_Space_Handler

  /*!
    \class Phase_Space_Handler
    \brief the main steering class for integration and event generation
  */

  /*!
    \var ATOOLS::Integration_Info *Phase_Space_Handler::p_info

    Phase_Space_Handler distributes information on the various 
    integration variables and weights through an instance of Integration_Info.
    Each Single_Channel is able to gain access to the variables via the assignment 
    of an Info_Key during its initialization.
  */

}//end of namespace PHASIC

#endif
