#ifndef SINIC_Beam_Remnants_Singlet_Sorter_H
#define SINIC_Beam_Remnants_Singlet_Sorter_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>
#include <map>

namespace SINIC {
  class Rap_Order {
  public :
    inline bool operator()(ATOOLS::Particle * p1,ATOOLS::Particle * p2)
    {
      if (p1->Momentum().Y()<p2->Momentum().Y()) return true;
      return false;
    }
  };
  typedef std::list<ATOOLS::Particle *> PartList;
  typedef std::set<ATOOLS::Particle *,Rap_Order> PartYSet;
  class Rap_Order2 {
  public :
    inline bool operator()(PartYSet * ps1,PartYSet * ps2) {
      if (ps1->empty()) return true;
      if (ps2->empty()) return false;
      if ((*ps1->begin())->Momentum().Y()<
	  (*ps2->begin())->Momentum().Y()) return true;
      return false;
    }
  };
  typedef std::set<PartYSet *,Rap_Order2> PSetYSet;
  typedef std::pair<ATOOLS::Particle *,
		    std::pair<unsigned int, unsigned int> > PCPair;
  typedef std::list<PCPair> PCList;

  class Singlet_Sorter {
  private:
    PartList * m_inlist;
    ATOOLS::Particle * FindNextStart();
    ATOOLS::Particle * FindNext(const unsigned int & ref1,
				const unsigned int & ref2=-1);
  public:
    Singlet_Sorter();
    ~Singlet_Sorter() {}
    void Sort(PartList * inlist,PSetYSet * outlists);
    void Sort(PartList * inlist,PCList * outlist);
  };
}
#endif

