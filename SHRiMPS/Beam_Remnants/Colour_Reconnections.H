#ifndef SHRIMPS_Beam_Remnants_Colour_Reconnections_H
#define SHRIMPS_Beam_Remnants_Colour_Reconnections_H

//#include "SHRiMPS/Beam_Remnants/Singlet_Sorter.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Return_Value.H"
#include <map>

namespace SHRIMPS {

  struct partcomp{
    bool operator() (const ATOOLS::Particle * part1, const ATOOLS::Particle * part2) const
    {
      if (part1->Number() < part2->Number()) return true;
      else return false;
    }
  };

  struct blobcomp{
    bool operator() (const ATOOLS::Blob * blob1, const ATOOLS::Blob * blob2) const
    {
      if (blob1->Id() < blob2->Id()) return true;
      else return false;
    }
  };

  typedef std::pair<unsigned int,unsigned int> colpair;
  typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> partpair;
  typedef std::map<double,ATOOLS::Particle *> partdists; 
  class Colour_Reconnections {
  private:
    bool           m_on, m_shuffled;
    double         m_reconn, m_Q02, m_b02, m_eta, m_smin, m_ycut;
    
    std::set<int>                m_colours;
    std::set<ATOOLS::Particle *, partcomp> m_trips, m_antis;
    std::map<ATOOLS::Particle *, partdists, partcomp> m_links;
    std::map<ATOOLS::Particle *, colpair, partcomp> m_newcols;
    std::list<partpair> m_pairs;

    bool                                        m_analyse;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    void HarvestParticles(ATOOLS::Blob_List * blobs);
    void FillWeightTable();
    void OutputWeightTable();
    void ShuffleColours();
    void SaveLastGluon(ATOOLS::Particle * part);
    double Distance(ATOOLS::Particle * part1,ATOOLS::Particle * part2,
			  const bool & spat=true);

    ATOOLS::Blob * AddReconnectionBlob();
    void AddParticleToReconnectionBlob(ATOOLS::Blob * blob,
				       ATOOLS::Particle * part);

    size_t ColourConnected(ATOOLS::Particle * part1,ATOOLS::Particle * part2);
  public:
    Colour_Reconnections();
    ~Colour_Reconnections();
    bool FinishConfiguration(ATOOLS::Blob_List * blobs,
			     const double & smin=-1.);
  };

  inline size_t Colour_Reconnections::
  ColourConnected(ATOOLS::Particle * part1,ATOOLS::Particle * part2) {
    for (size_t i=1;i<3;i++) {
      if (part1->GetFlow(i)==part2->GetFlow(3-i)) return i;
    }
    return 0;
  }
}
#endif

