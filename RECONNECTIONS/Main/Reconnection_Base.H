#ifndef RECONNECTIONS_Main_Reconnection_Base_H
#define RECONNECTIONS_Main_Reconnection_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Return_Value.H"

namespace RECONNECTIONS {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Reconnection_Base {
  protected:
    ATOOLS::ParticleSet                        m_parts[2];
    std::map<unsigned int, ATOOLS::Particle *> m_cols[2];
    ATOOLS::Part_List                          m_particles;
    bool                                       m_found;
    
    bool                                        m_analysis;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    virtual void SetParameters() = 0;
    bool HarvestParticles(ATOOLS::Blob_List * blobs);
    void HarvestParticleInfo(ATOOLS::Particle * part);
    void FillMassesInHistogram(ATOOLS::Histogram * histo);

    std::vector<double> variation_weights;
    std::vector<double> tmp_variation_weights;

  public:
    Reconnection_Base();
    virtual ~Reconnection_Base();

    virtual void   Initialize();
    virtual void   Reset();
    virtual int    operator()(ATOOLS::Blob_List *const blobs) = 0;
    virtual double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti) = 0;

    ATOOLS::Part_List              * GetParticles()               { return &m_particles; }
    ATOOLS::ParticleSet            * GetParts(const size_t & pos) { return &m_parts[pos]; }

    void reset_var_weights() {
      std::fill(tmp_variation_weights.begin(), tmp_variation_weights.end(), 1);
    }

    std::vector<double>& get_variationweights() {
      return variation_weights;
    };

    void reset_variationweights(int nvars) {
      //std::cout << "resetting variation weights" << std::endl;
      //std::cout << "---------------------------" << std::endl;
      tmp_variation_weights.resize(nvars);
      reset_var_weights();
      
      variation_weights.resize(nvars);
      std::fill(variation_weights.begin(), variation_weights.end(), 1);
    };

  };
}
#endif
