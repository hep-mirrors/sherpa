#ifndef PHASIC_Channels_RS_Multi_Channel_H
#define PHASIC_Channels_RS_Multi_Channel_H

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/RS_Emission_Generator.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace PHASIC {

  class Process_Base;
  class Phase_Space_Handler;

  class RS_Multi_Channel: public Multi_Channel {
  private:

    Process_Base *p_proc;

    Multi_Channel *p_fsmc;

    RS_Emission_Generator m_eeg;

  public:

    RS_Multi_Channel(Process_Base *const proc,
		     Phase_Space_Handler *const psh);

    ~RS_Multi_Channel();

    void Reset();
    bool Initialize();
    void Add(Single_Channel *channel);  

    void Optimize(double);
    void EndOptimize(double);
    void MPISync();

    void AddPoint(double value);

    void GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts,
			bool compute=false);
    void GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts);
    
    void WriteOut(std::string pId);
    bool ReadIn(std::string pId);

    size_t Number();

    std::string Name();     
    std::string ChID();

    bool OptimizationFinished();

    void ISRInfo(int i,int &t,double &m,double &w);
    void ISRInfo(std::vector<int> &ts,std::vector<double> &ms,
		 std::vector<double> &ws) const;

    void Print();

  };// end of class RS_Multi_Channel

}// end of namespace PHASIC

#endif
