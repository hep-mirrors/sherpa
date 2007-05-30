#ifndef PreSample_Multi_Channel_h
#define PreSample_Multi_Channel_h

#include "Multi_Channel.H"
#include "Vegas.H"
#include "Color_Integrator.H"

namespace PHASIC {

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
    PreSample_Multi_Channel(std::string,Color_Integrator *const colint);
    ~PreSample_Multi_Channel();
    void AddMC(Multi_Channel *);  
    
    void DropAllChannels();     
    std::vector<Single_Channel*> ExtractChannels();
    void Reset();
    void ResetOpt();
    void Optimize(double);
    void EndOptimize(double);
    bool OptimizationFinished(); 
    void AddPoint(double);

    void GeneratePoint(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts);
    void GenerateWeight(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts); 
 
    double GetAlpha(int type)       { return m_alpha[type]; }
    size_t              Number()    { return m_subchannels.size(); }
    void Print();
    void WriteOut(std::string);

    int                 HandicapFactor();
  };



}


#endif
