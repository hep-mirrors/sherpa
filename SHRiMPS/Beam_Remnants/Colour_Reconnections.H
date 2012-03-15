#ifndef SHRIMPS_Beam_Remnants_Colour_Reconnections_H
#define SHRIMPS_Beam_Remnants_Colour_Reconnections_H

#include "SHRiMPS/Beam_Remnants/Singlet_Sorter.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Return_Value.H"
#include <map>

namespace SHRIMPS {
  class Colour_Reconnections {
  private:
    bool           m_on, m_shuffled;
    double         m_colfac, m_reconn, m_Q02, m_b02, m_inveta, m_smin, m_ycut;
    PCList         m_pclist;
    Singlet_Sorter m_sorter;
    
    bool                                        m_analyse;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    void HarvestParticles(ATOOLS::Blob_List * blobs);
    void ShuffleColours();
    void SkewList(PCList::iterator & pit1,PCList::iterator & pit2,
		  PCList::iterator & pit3,PCList::iterator & pit4);
    const double Weight(ATOOLS::Particle * part1,ATOOLS::Particle * part2,
			const bool & spatial=true);

    ATOOLS::Blob * AddReconnectionBlob();

    bool ColourConnected(ATOOLS::Particle * part1,ATOOLS::Particle * part2);
  public:
    Colour_Reconnections();
    ~Colour_Reconnections();
    bool FinishConfiguration(ATOOLS::Blob_List * blobs,
			     const double & smin=-1.);
  };

  inline bool Colour_Reconnections::
  ColourConnected(ATOOLS::Particle * part1,ATOOLS::Particle * part2) {
    for (size_t i=1;i<3;i++) {
      if (part1->GetFlow(i)==part2->GetFlow(3-i)) return true;
    }
    return false;
  }
}
#endif

