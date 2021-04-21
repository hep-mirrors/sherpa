#ifndef RECONNECTIONS_Main_Reconnect_By_Singlet_H
#define RECONNECTIONS_Main_Reconnect_By_Singlet_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"
#include "RECONNECTIONS/Main/Reconnection_Weights.H"

namespace RECONNECTIONS {
  class Reconnect_By_Singlet : public Reconnection_Base {
  private:    
    Reconnection_Weights m_weights;
    int    m_Pmode;
    double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_restring;
    std::list<ATOOLS::Part_List *> m_singlets;    
    
    void SetParameters();
    void MakeSinglets();
    void FlattenSinglets();
    ATOOLS::Particle * FindStart();
    ATOOLS::Particle * FindNext(const size_t & col);
    void ReshuffleSinglets();
    bool ReshuffleSinglet(ATOOLS::Part_List * singlet);
    void ReconnectSinglets();
    void SpliceSinglets(ATOOLS::Part_List * sing1,ATOOLS::Part_List * sing2,
			ATOOLS::Part_Iterator & pit1,ATOOLS::Part_Iterator & pit2);
    void AftermathOfSlicing(ATOOLS::Particle * part11,ATOOLS::Particle * part12,
			    ATOOLS::Particle * part21,ATOOLS::Particle * part22);
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  public:
    Reconnect_By_Singlet();
    ~Reconnect_By_Singlet();

    bool   operator()(ATOOLS::Blob_List *const blobs);
    void   Reset();
    double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  };
}
#endif
