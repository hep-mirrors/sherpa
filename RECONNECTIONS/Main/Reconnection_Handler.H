#ifndef RECONNECTIONS_Main_Reconnection_Handler_H
#define RECONNECTIONS_Main_Reconnection_Handler_H

#include "RECONNECTIONS/Main/Reconnection_Weights.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Return_Value.H"

namespace RECONNECTIONS {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Reconnection_Handler {
  private:
    bool m_on;
    ATOOLS::ParticleSet                        m_parts[2];
    std::map<unsigned int, ATOOLS::Particle *> m_cols[2];
    Reconnection_Weights                       m_weights;
    std::list<ATOOLS::Part_List *>             m_singlets;

    bool                                        m_analysis;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    bool HarvestParticles(ATOOLS::Blob_List * blobs);
    void HarvestParticleInfo(ATOOLS::Particle * part);
    void MakeSinglets();
    ATOOLS::Particle * FindStart();
    ATOOLS::Particle * FindNext(const size_t & col);
    void ReshuffleSinglets();
    bool ReshuffleSinglet(ATOOLS::Part_List * singlet);
    void ReconnectSinglets();
    void SpliceSinglets(ATOOLS::Part_List * sing1,ATOOLS::Part_List * sing2,
			ATOOLS::Part_Iterator & pit1,ATOOLS::Part_Iterator & pit2);
    void AftermathOfSlicing(ATOOLS::Particle * part11,ATOOLS::Particle * part12,
			    ATOOLS::Particle * part21,ATOOLS::Particle * part22);
    void AddReconnectionBlob(ATOOLS::Blob_List * blobs);
    void FillMassesInHistogram(ATOOLS::Histogram * histo);
 public:
    Reconnection_Handler();
    ~Reconnection_Handler();
    
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);

    void Initialize(ATOOLS::Default_Reader *const reader);
    void Reset();

    ATOOLS::ParticleSet * const GetParts(const size_t & pos) { return &m_parts[pos]; }
    std::list<ATOOLS::Part_List *> * const GetSinglets() { return &m_singlets; }
  };
}

#endif
