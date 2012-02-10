#ifndef SINIC_Beam_Remnants_Colour_Generator_H
#define SINIC_Beam_Remnants_Colour_Generator_H

#include "SINIC++/Beam_Remnants/Hadron_Dissociation.H"
#include "SINIC++/Event_Generation/Ladder.H"
#include "ATOOLS/Phys/Blob.H"

namespace SINIC {
  class Colour_Generator {
  private:
    std::vector<Hadron_Dissociation * > m_hadrons;

    std::set<int>      m_col[2][2];
    ATOOLS::Particle * p_inpart1, * p_inpart2;
    ATOOLS::Blob     * p_blob, * p_compensator;
    Ladder           * p_ladder;
    
    bool   Primary(ATOOLS::Particle ** parts,const size_t & counter);
    bool   Rescatter();
    bool   FixUncorrelatedIndices();
    bool   FixCorrelatedIndices();
    size_t FixAColourLine();
    size_t SelectColourReplacement(ATOOLS::Particle * part1,
				   ATOOLS::Particle * part2);
    int    ModifyOriginators(ATOOLS::Particle * part1,
			     ATOOLS::Particle * part2);
    bool   ReplaceColour(ATOOLS::Blob * blob,const size_t & index,
			 const size_t & col,const size_t & col1,
			 const size_t & adj);
    int    CorrelatedIndices(int col[2][2]);
    int    UncorrelatedIndices(int col[2][2]);
    int    PickOneColour(const size_t & beam,const size_t & index,
			 const int & avoid=-1);
    int    PickColourPair(const size_t & beam,const size_t & index);
    size_t PickIndexAndColour(int & col);
    void   PickTwoColours(const size_t & beam,int * cols);
  public:
    Colour_Generator(std::vector<Hadron_Dissociation * > * hadrons);
    ~Colour_Generator();

    bool operator()(Ladder * ladder,ATOOLS::Particle ** parts,
		    const size_t & counter);
    void FinalColours();

    inline void Reset() { 
      for (size_t beam=0;beam<2;beam++) { 
	for (size_t flow=0;flow<2;flow++) { 
	  m_col[beam][flow].clear();
	} 
      } 
    }
    inline void SetSoftBlob(ATOOLS::Blob * blob) { p_blob = blob; }
    inline void SetInitials(ATOOLS::Particle * part1,ATOOLS::Particle * part2) {
      p_inpart1 = part1; p_inpart2 = part2;
    }
    inline ATOOLS::Blob * GetCompensatorBlob() const { return p_compensator; }
    inline void PrintPendingColours() {
      for (size_t beam=0;beam<2;beam++) {
	msg_Out()<<"  Cols in beam("<<beam<<"): { ";
	for (size_t index=0;index<2;index++) {
	  if (index==1) msg_Out()<<" { ";
	  for (std::set<int>::iterator cit=m_col[beam][index].begin();
	       cit!=m_col[beam][index].end();cit++) 
	    msg_Out()<<(*cit)<<" ";
	  msg_Out()<<"}";
	}
	msg_Out()<<"\n";
	m_hadrons[1-beam]->PrintParticles();
      }
    }
    inline size_t ColourConnected(ATOOLS::Particle * part1,
				  ATOOLS::Particle * part2) const {
      size_t index(0);
      for (size_t i=1;i<3;i++) {
	if (part1->GetFlow(i)==part2->GetFlow(3-i) &&
	    part1->GetFlow(i)!=0) index+=i;
      }
      return index;
    }
  };
}

#endif
