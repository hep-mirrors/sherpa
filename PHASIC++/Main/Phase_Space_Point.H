#ifndef PHASIC_Main_Phase_Space_Point_H
#define PHASIC_Main_Phase_Space_Point_H

#include "PHASIC++/Channels/Beam_Channels.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <iostream>

namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF  { class ISR_Handler; }

namespace PHASIC {
  class Phase_Space_Handler;
  class Process_Integrator;
  class Cut_Data;

  struct psmode {    
    enum code {
      normal_call = 0,
      no_lim_isr  = 32,
      no_gen_isr  = 64
    };
  };// end of struct psm
  inline psmode::code operator|(const psmode::code &c1,const psmode::code &c2) {
    return (psmode::code)((int)c1|(int)c2);
  }
  inline psmode::code operator&(const psmode::code &c1,const psmode::code &c2) {
    return (psmode::code)((int)c1&(int)c2);
  }

  
  class Phase_Space_Point {
  private:
    Phase_Space_Handler             * p_pshandler;
    BEAM::Beam_Spectra_Handler      * p_beamhandler;
    PDF::ISR_Handler                * p_isrhandler;
    ATOOLS::Vec4D                   * p_moms;
    Cut_Data                        * p_cuts;
    Multi_Channel  * p_beamchannels, * p_isrchannels, * p_fsrchannels;
    ATOOLS::Info_Key m_beamspkey, m_beamykey, m_isrspkey, m_isrykey;
    
    psmode::code        m_mode;    
    size_t              m_nin, m_nout, m_nvec;
    ATOOLS::Poincare  * p_fixedISboost;
    double              m_Ecms, m_Eprime, m_smin, m_sprime, m_y, m_fixedsprime, m_fixedy;
    std::vector<double> m_masses;
    double              m_osmass, m_masses2[2];
    ATOOLS::Vec4D       m_ISmoms[2];
    double              m_weight, m_ISsymmetryfactor;

    void   InitFixedIncomings();
    bool   DefineBeamKinematics();
    bool   DefineISRKinematics(Process_Integrator *const process);
    bool   DefineFSRKinematics();
    void   BoostSystem();
    bool   Check4Momentum();
    void   CorrectMomenta();

    inline void Reset(const psmode::code & mode) {
      m_mode   = mode;
      m_sprime = m_fixedsprime; m_y = m_fixedy;
      m_ISsymmetryfactor = 1.;
    }
  public:
    Phase_Space_Point();
    ~Phase_Space_Point();

    void   Init(Phase_Space_Handler * psh);
    void   InitCuts(Process_Integrator *const process);
    bool   operator()(Process_Integrator *const process,
		      const psmode::code & mode=psmode::normal_call);

    double CalculateWeight();
    bool   MakeIncoming();
    
    inline void SetBeamIntegrator(Multi_Channel * channels) { p_beamchannels = channels; }
    inline void SetISRIntegrator(Multi_Channel * channels)  { p_isrchannels  = channels; }
    inline void SetFSRIntegrator(Multi_Channel * channels)  { p_fsrchannels  = channels; }
    inline void SetMomenta(ATOOLS::Vec4D * moms) {
      p_moms = moms;
      ATOOLS::Vec4D pin = p_moms[0]+p_moms[1];
      m_isrspkey[3] = pin.Abs2();
      m_isrykey[2]  = pin.Y();
    }
    inline void SetOSMass(const double &osmass)             { m_osmass       = osmass; }

    inline Multi_Channel * BeamIntegrator()   const { return p_beamchannels; }
    inline Multi_Channel * ISRIntegrator()    const { return p_isrchannels; }
    inline Multi_Channel * FSRIntegrator()    const { return p_fsrchannels; }
    inline Cut_Data      * Cuts()             const { return p_cuts; }
    inline const double  & Weight()           const { return m_weight; }
    inline const double  & ISSymmetryFactor() const { return m_ISsymmetryfactor; }


    inline void AddPoint(const double & value) {
      if (p_beamchannels) p_beamchannels->AddPoint(value);
      if (p_isrchannels)  p_isrchannels->AddPoint(value);
      p_fsrchannels->AddPoint(value);
    }
    inline void Optimize(const double & error) {
      if (p_beamchannels) p_beamchannels->Optimize(error);
      if (p_isrchannels)  p_isrchannels->Optimize(error);
      p_fsrchannels->Optimize(error);
    }
    inline void EndOptimize(const double & error) {
      if (p_beamchannels) p_beamchannels->EndOptimize(error);
      if (p_isrchannels)  p_isrchannels->EndOptimize(error);
      p_fsrchannels->EndOptimize(error);
    }

    inline void MPISync() {
      if (p_beamchannels) p_beamchannels->MPISync();
      if (p_isrchannels) p_isrchannels->MPISync();
      p_fsrchannels->MPISync();
    }
    inline void WriteOut(const std::string &pID) {
      if (p_beamchannels) p_beamchannels->WriteOut(pID+"/MC_Beam");
      if (p_isrchannels) p_isrchannels->WriteOut(pID+"/MC_ISR");
      if (p_fsrchannels) p_fsrchannels->WriteOut(pID+"/MC_FSR");
    }
    inline bool ReadIn(const std::string &pID,const size_t exclude) {
      bool okay = true;
      if (p_beamchannels && !(exclude&1)) okay = okay && p_beamchannels->ReadIn(pID+"/MC_Beam");
      if (p_isrchannels && !(exclude&2))  okay = okay && p_isrchannels->ReadIn(pID+"/MC_ISR");
      if (p_fsrchannels && !(exclude&16)) okay = okay && p_fsrchannels->ReadIn(pID+"/MC_FSR");
      return okay;
    }    
    
    
    void Print(std::ostream & str);
  };
}
#endif
