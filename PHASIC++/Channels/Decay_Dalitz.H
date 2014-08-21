#ifndef PHASIC_Channels_Decay_Dalitz_H
#define PHASIC_Channels_Decay_Dalitz_H

#include "PHASIC++/Channels/Single_Channel.H"

namespace PHASIC {
  class Decay_Dalitz : public Single_Channel {
    
  private:
    ATOOLS::Vec4D m_decvec, m_pvec;
    double        m_pmass, m_pwidth, m_sexp, m_smin, m_smax;
    int           m_p1, m_p2, m_dir, m_mode;
    std::string   m_propmass;				
    int           m_chnumber;			
    const ATOOLS::Mass_Selector* p_masssel;
  public :
    Decay_Dalitz(const ATOOLS::Flavour * fl,
                 const double& mass, const double& width,
                 size_t dir, size_t p1, size_t p2,
                 const ATOOLS::Mass_Selector* masssel);
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);  
    std::string Name() { return std::string("Decay_Dalitz"); }
    void MPISync() {}
  };
}

#endif
