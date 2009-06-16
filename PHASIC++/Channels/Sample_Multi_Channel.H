#ifndef PHASIC_Channels_Sample_Multi_Channel_H
#define PHASIC_Channels_Sample_Multi_Channel_H

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Channels/Vegas.H"

#include <map>

namespace PHASIC {

  class Color_Integrator;

  typedef std::vector<Single_Channel*> Channel_Vector;

  struct Idx_Key {
    size_t     m_id;
    Idx_Vector m_t, m_b;
    Idx_Key(): m_id(0) {}
    Idx_Key(const size_t &id,const Idx_Vector &t,
	    const Idx_Vector &b): m_id(id), m_t(t), m_b(b) {}
  };// end of struct Idx_Key

  std::ostream &operator<<(std::ostream &str,const Idx_Key &key);

  typedef std::map<Idx_Vector,Idx_Key> Index_Map;
  typedef std::vector<Idx_Key>         Key_Vector;

  typedef std::vector<Single_Channel*> Channel_Vector;

  class Phase_Space_Handler;

  class Sample_Multi_Channel: public Multi_Channel {
  private:

    Color_Integrator *p_colint;

    Channel_Vector m_channels, m_oldchannels;
    Key_Vector     m_confs;
    Double_Vector  m_alpha;
    Index_Map      m_idx;
    ATOOLS::Vec4D_Vector   m_mp;

    size_t  m_cur, m_first;

    Vegas  *p_vegas;
    double *p_ran, m_wsum;

    bool m_otfcc;

    Idx_Key GenerateKey(Idx_Vector &perm);

    double SingleWeight
    (ATOOLS::Vec4D *p,Cut_Data *cuts,const Idx_Key &key);

    void PDWeight(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void PDPoint(ATOOLS::Vec4D *p,Cut_Data *cuts);

    void OTFWeight(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void OTFPoint(ATOOLS::Vec4D *p,Cut_Data *cuts);

    void FillIdxMap();

  public:

    Sample_Multi_Channel(Phase_Space_Handler *const psh,
			 Color_Integrator *const ci);

    ~Sample_Multi_Channel();

    bool Initialize(const Channel_Vector &chs);
    
    void Clear();
    void Reset();

    void Optimize(double);
    void EndOptimize(double);

    void AddPoint(double value);

    void GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts);
    void GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts);
 
    void WriteOut(std::string pId);
    bool ReadIn(std::string pId);

    std::string Name();     
    std::string ChID();

    void ISRInfo(int &,double &,double &);

    void Print();

    inline void SetValidN(const long unsigned int &n) { n_contrib=n; }
    inline void SetN(const long unsigned int &n)      { n_points=n;  }

  };// end of class Sample_Multi_Channel

}// end of namespace PHASIC

#endif
