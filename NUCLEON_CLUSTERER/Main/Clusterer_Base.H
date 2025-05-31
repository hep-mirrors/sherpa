#ifndef NUCLEON_CLUSTERER_Main_Clusterer_Base_H
#define NUCLEON_CLUSTERER_Main_Clusterer_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Return_Value.H"

namespace NUCLEON_CLUSTERER {
  typedef std::pair<unsigned int,unsigned int> nucleonpair;
  
  class Clusterer_Base {
  protected:
    ATOOLS::ParticleSet                        m_parts[2]; // TODO: change this to show nucleon pairs or anti nucleon pairs
    std::map<unsigned int, ATOOLS::Particle *> m_cols[2]; // TODO: change away from colour stuff
    ATOOLS::Part_List                          m_particles; // all particles in the cluster
    bool                                       m_found; // true if particles were found
    
    bool                                        m_analysis;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    virtual void SetParameters() = 0;
    bool HarvestParticles(ATOOLS::Blob_List * blobs); // collect particles from event blobs
    void HarvestParticleInfo(ATOOLS::Particle * part); // process inidivual particles/ in the future nucleons
    // bool BalanceColours(); // here we would check for colour but for nucleons check isospin TODO
    void FillMassesInHistogram(ATOOLS::Histogram * histo);
  public:
    Clusterer_Base();
    virtual ~Clusterer_Base();

    virtual void   Initialize();
    virtual void   Reset();
    // main clustering logic here TODO 
    virtual int    operator()(ATOOLS::Blob_List *const blobs) = 0;
    // define clustering metric here TODO 
    virtual double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) = 0;

    ATOOLS::Part_List              * GetParticles()               { return &m_particles; }
    ATOOLS::ParticleSet            * GetParts(const size_t & pos) { return &m_parts[pos]; }
  };
}
#endif
