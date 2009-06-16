#ifndef PHASIC_Channels_PreSample_Multi_Channel_H
#define PHASIC_Channels_PreSample_Multi_Channel_H

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Main/Color_Integrator.H"

namespace PHASIC {

  class Phase_Space_Handler;

  class PreSample_Multi_Channel : public Multi_Channel {
    std::vector<Multi_Channel *> m_subchannels;
    std::vector<double> m_alpha;
    std::vector<double> m_palpha;
    Vegas* p_vegas;
    double* p_ran;
    double m_calpha;
    std::vector<double> m_multi;
    Color_Integrator *p_colint;
  public:
    PreSample_Multi_Channel(Phase_Space_Handler *const psh,
			    Color_Integrator *const colint);
    ~PreSample_Multi_Channel();
    void AddMC(Multi_Channel *);  
    
    void DropAllChannels(const bool del=true);     
    std::vector<Single_Channel*> ExtractChannels();
    void Reset();
    void ResetOpt();
    void Optimize(double);
    void EndOptimize(double);
    bool OptimizationFinished(); 
    void AddPoint(double);

    void GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts); 
 
    double GetAlpha(int type)       { return m_alpha[type]; }
    size_t              Number()    { return m_subchannels.size(); }
    void Print();
    void WriteOut(std::string);
    bool ReadIn(std::string);

    int                 HandicapFactor();
  };



}


#endif
