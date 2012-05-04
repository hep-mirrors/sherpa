#ifndef SHRIMPS_Beam_Remnants_Rescatter_Handler_H
#define SHRIMPS_Beam_Remnants_Rescatter_Handler_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Beam_Remnants/Singlet_Sorter.H"
#include "MODEL/Main/Strong_Coupling.H"

namespace SHRIMPS {
  typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> PartPair;
  typedef std::map<double,PartPair,std::greater<double> >  ProbPairMap;

  class Rescatter_Handler {
  private:
    bool                      m_rescatter, m_mustmatch;
    double                    m_rescprob, m_singprob;
    Beam_Remnant_Handler    * p_beams;
    MODEL::Strong_Coupling  * p_alphaS;
    Singlet_Sorter            m_sorter;
    double                    m_Ylimit;

    bool                                        m_analyse;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    Omega_ik                            * p_eikonal;
    int                                   m_Nresc, m_Nfact;
    std::set<ATOOLS::Blob * >             m_treatedblobs;
    std::set<ATOOLS::Particle * >         m_particles;
    std::list<std::pair<double, double> > m_intervals;
    ProbPairMap                           m_probpairs;
    double                                m_smin,m_B,m_b1,m_b2;
 

    std::map<ATOOLS::Particle *,ATOOLS::Particle *> m_partmap;

    bool DealWithBlob(ATOOLS::Blob * blob);
    void AddParticleToRescatters(ATOOLS::Particle * part);
    void DeleteProbPairs(ATOOLS::Particle *& lpart1,ATOOLS::Particle *& lpart2);
    bool CanRescatter(const ATOOLS::Particle * lpart1,
		      const ATOOLS::Particle * lpart2) const;
    bool IsColourConnected(const ATOOLS::Particle * lpart1,
			   const ATOOLS::Particle * lpart2) const;
  public:
    Rescatter_Handler() {}
    Rescatter_Handler(Beam_Remnant_Handler * beams);
    ~Rescatter_Handler();

    void   ResetCollision(Omega_ik * eikonal,const double & smin,
			  const double & B);
    void   UpdateCollision(ATOOLS::Blob_List * blobs);
    void   ResetRescatter(const bool & enforce=true);
    bool   SelectRescatter(ATOOLS::Particle *& part1,ATOOLS::Particle *& part2);
    bool   ConnectBlobs(ATOOLS::Blob_List * blobs,ATOOLS::Blob * add);
    void   FillInitialStateIntoBlob(ATOOLS::Blob * blob,Ladder * ladder);

    inline void Map(ATOOLS::Particle * p1,ATOOLS::Particle * p2) {
      if (p1->Number()<p2->Number()) m_partmap[p1] = p2; 
      else m_partmap[p2] = p1;      
    }
    inline ATOOLS::Particle * FindMap(ATOOLS::Particle * part) {
      return (m_partmap.find(part)!=m_partmap.end()?m_partmap[part]:NULL);
    }
  };
}
#endif
