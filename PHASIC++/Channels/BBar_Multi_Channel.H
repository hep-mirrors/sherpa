#ifndef PHASIC_Channels_BBar_Multi_Channel_H
#define PHASIC_Channels_BBar_Multi_Channel_H

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/BBar_Emission_Generator.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace PHASIC {

  class Process_Base;
  class Phase_Space_Handler;

  class BBar_Multi_Channel: public Multi_Channel {
  private:

    Process_Base *p_proc;

    Multi_Channel *p_fsmc;
    Cut_Data      *p_cuts;

    BBar_Emission_Generator m_eeg;

  public:

    BBar_Multi_Channel(Process_Base *const proc,
		       Process_Base *const sproc,
		       Phase_Space_Handler *const psh);

    ~BBar_Multi_Channel();

    Dipole_Params Active(Process_Base *const bviproc) const;

    void Reset();

    void Optimize(double);
    void EndOptimize(double);
    void MPISync();

    void AddPoint(double value);

    void GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts);
    
    void GenerateEmissionPoint
    (const ATOOLS::Cluster_Amplitude &ampl);

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

  };// end of class BBar_Multi_Channel

}// end of namespace PHASIC

#endif
