#ifndef PHASIC_Channels_POWHEG_Multi_Channel_H
#define PHASIC_Channels_POWHEG_Multi_Channel_H

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Extra_Emission_Generator.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace PHASIC {

  class Process_Base;
  class POWHEG_Process;
  class Phase_Space_Handler;

  class POWHEG_Multi_Channel: public Multi_Channel {
  private:

    POWHEG_Process *p_proc;
    Multi_Channel  *p_fsmc, *p_ismc;

    Extra_Emission_Generator m_eeg;

    ATOOLS::Vec4D_Vector       m_pb, m_pr;
    ATOOLS::Cluster_Amplitude *p_bampl;

    Cut_Data         *p_cuts;
    ATOOLS::Info_Key  m_isrspkey, m_isrykey;

    double m_bweight;
    size_t m_ismode, m_emode;

  public:

    POWHEG_Multi_Channel(POWHEG_Process *const proc,
			 Process_Base *const sproc,
			 Phase_Space_Handler *const psh);

    ~POWHEG_Multi_Channel();

    void Reset();

    void Optimize(double);
    void EndOptimize(double);
    void MPISync();

    void AddPoint(double value);

    void GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts);
    
    void GenerateEmissionPoint(const ATOOLS::Vec4D_Vector& pb);
 
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

    double GenerateWeight(const ATOOLS::Cluster_Amplitude &ampl);

    size_t ActiveIdx() const;

    double SelectionWeight(const size_t &idx) const;

    inline double BRWeight() const { return m_bweight/m_weight; }

    inline ATOOLS::Cluster_Amplitude *BAmpl() const { return p_bampl; }
    inline const ATOOLS::Vec4D_Vector& RealMoms() const { return m_pr; }

  };// end of class POWHEG_Multi_Channel

}// end of namespace PHASIC

#endif
